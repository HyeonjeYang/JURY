"""Audit a sequencing primer walk and design primers for the gaps."""
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt, gc_fraction
from .seqio import revcomp


class Primer:
    __slots__ = ("name", "start", "end", "seq", "strand")

    def __init__(self, name, start, end, seq, strand=1):
        self.name, self.start, self.end = name, start, end
        self.seq, self.strand = seq, strand

    @property
    def tm(self):
        return mt.Tm_NN(Seq(self.seq))

    @property
    def gc(self):
        return 100 * gc_fraction(self.seq)

    @property
    def three_prime(self):
        """Map coordinate the read starts extending from."""
        return self.end if self.strand >= 0 else self.start

    def __repr__(self):
        return f"<{self.name} {self.start}-{self.end} {len(self.seq)}nt Tm{self.tm:.1f}>"


def occurrences(seq, pmap):
    """(forward, reverse) hit counts on the circular map."""
    n, d = len(pmap), pmap.doubled
    s = str(seq).upper()
    if not s or len(s) > n:
        return 0, 0
    search = d[:n + len(s) - 1]

    def count(motif):
        hits, offset = 0, 0
        while True:
            offset = search.find(motif, offset)
            if offset < 0:
                return hits
            hits += 1
            offset += 1

    return count(s), count(revcomp(s))


def has_run(seq, n=4):
    c, prev = 1, ""
    for b in seq:
        c = c + 1 if b == prev else 1
        prev = b
        if c >= n:
            return True
    return False


def self_dimer3(seq, n=5):
    return revcomp(seq[-n:]) in seq


def collect_primers(pmap, pattern=r"seq_primer"):
    """Pull annotated sequencing primers off the map."""
    out = []
    for f in pmap.features_named(pattern):
        seq = pmap.slice(f.start, f.end)
        if f.strand < 0:
            seq = revcomp(seq)
        out.append(Primer(f.name, f.start, f.end, seq, f.strand))
    out.sort(key=lambda p: p.three_prime)
    return out


def audit(primers, pmap, dead=50, good=700):
    """Model each read window and report plasmid-wide coverage and gaps.

    dead: bases after the 3' end that are not yet reliable.
    good: last reliable base, counted from the 3' end.
    """
    n = len(pmap)
    if n == 0:
        raise ValueError("cannot audit an empty plasmid map")
    if dead < 0 or good < dead:
        raise ValueError("read window must satisfy 0 <= dead <= good")
    cov = [0] * (n + 1)
    windows = []
    for p in primers:
        direction = 1 if p.strand >= 0 else -1
        a = p.three_prime + direction * dead
        b = p.three_prime + direction * good
        windows.append((p, (a - 1) % n + 1, (b - 1) % n + 1))
        for offset in range(dead, good + 1):
            q = (p.three_prime + direction * offset - 1) % n + 1
            cov[q] += 1
    gaps, run = [], None
    for q in range(1, n + 1):
        if cov[q] == 0 and run is None:
            run = q
        elif cov[q] > 0 and run is not None:
            gaps.append((run, q - 1)); run = None
    if run is not None:
        gaps.append((run, n))
    missing = sum(1 for q in range(1, n + 1) if cov[q] == 0)
    return {"windows": windows, "gaps": gaps, "uncovered": missing,
            "pct": 100.0 * (n - missing) / n, "depth": cov}


def score_candidate(seq, pmap, tm_lo=54, tm_hi=63, gc_lo=40, gc_hi=62):
    """None if the candidate fails a hard filter, else a preference score."""
    if len(seq) < 2 or set(seq) - set("ACGT"):
        return None
    tm, gc = mt.Tm_NN(Seq(seq)), 100 * gc_fraction(seq)
    if not (tm_lo <= tm <= tm_hi) or not (gc_lo <= gc <= gc_hi):
        return None
    if has_run(seq) or self_dimer3(seq):
        return None
    f, r = occurrences(seq, pmap)
    if f != 1 or r != 0:
        return None
    s = -abs(tm - 58) * 2 - abs(gc - 50) * 0.3
    if seq[-1] in "GC":
        s += 2.0
    if seq[-2] in "GC":
        s += 0.5
    return s


def design_at(target, pmap, window=70, lengths=range(18, 25), accept=None, **kw):
    """Best forward primer with its 3' end near target."""
    n, best = len(pmap), None
    if window < 1:
        raise ValueError("primer search window must be positive")
    for end in range(target - window, target + window):
        for length in lengths:
            if length < 2 or length > n:
                continue
            start = (end - length) % n + 1
            finish = (end - 1) % n + 1
            seq = pmap.slice(start, start + length - 1)
            if len(seq) < length:
                continue
            candidate = Primer(f"new@{finish}", start, finish, seq, 1)
            if accept is not None and not accept(candidate):
                continue
            score = score_candidate(seq, pmap, **kw)
            if score is not None and (best is None or score > best[0]):
                best = score, candidate
    return None if best is None else best[1]


def fill_gaps(primers, pmap, spacing=600, progress=True, dead=50, good=700,
              **design_kw):
    """Design forward primers that improve the modeled sequencing coverage.

    Sequence filters can prevent closure; audit the returned walk for gaps.
    """
    if spacing <= 0:
        raise ValueError("primer spacing must be positive")
    primers = list(primers)
    n = len(pmap)
    state = audit(primers, pmap, dead, good)
    if not state["uncovered"]:
        return []
    ends = sorted({p.three_prime for p in primers}) or [1]
    targets = []
    for i, end in enumerate(ends):
        nxt = ends[(i + 1) % len(ends)]
        distance = (nxt - end) % n or n
        count = int(distance // spacing)
        for j in range(1, count + 1):
            targets.append((end + int(distance * j / (count + 1))) % n or n)
    bar = None
    if progress:
        try:
            from tqdm import tqdm
            bar = tqdm(desc="designing", unit="candidate")
        except ImportError:
            pass
    out = []
    seen = {(p.start, p.end, p.strand) for p in primers}

    def add(primer):
        nonlocal state
        if primer is None or (primer.start, primer.end, primer.strand) in seen:
            return False
        updated = audit(primers + out + [primer], pmap, dead, good)
        if updated["uncovered"] >= state["uncovered"]:
            return False
        seen.add((primer.start, primer.end, primer.strand))
        out.append(primer)
        state = updated
        return True

    try:
        for target in targets:
            if not state["uncovered"]:
                break
            add(design_at(target, pmap, **design_kw))
            if bar is not None:
                bar.update()

        stalled = set()
        while state["uncovered"]:
            anchors = [start for start, _ in state["gaps"] if start not in stalled]
            if not anchors:
                break
            anchor = anchors[0]

            def covers_anchor(primer):
                distance = (anchor - primer.three_prime) % n
                return any(dead <= distance + lap * n <= good
                           for lap in range(good // n + 1))

            primer = design_at(anchor - dead, pmap, accept=covers_anchor, **design_kw)
            if not add(primer):
                stalled.add(anchor)
            if bar is not None:
                bar.update()
        return out
    finally:
        if bar is not None:
            bar.close()
