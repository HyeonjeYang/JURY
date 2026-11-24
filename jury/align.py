"""Align Sanger reads to a circular plasmid map."""
from Bio import Align
from .seqio import revcomp


class Diff:
    """One difference, with coordinates in the oriented read and map."""
    __slots__ = ("ref", "read", "kind", "detail", "size")

    def __init__(self, ref, read, kind, detail, size=1):
        self.ref, self.read = int(ref), int(read)
        self.kind, self.detail, self.size = kind, detail, int(size)

    def __repr__(self):
        return f"<{self.kind} ref{self.ref} read{self.read} {self.detail}>"

    def positions(self, plasmid_len):
        """A deletion affects its entire interval, including across the origin."""
        size = self.size if self.kind == "del" else 1
        return {((self.ref - 1 + k) % plasmid_len) + 1 for k in range(size)}


class ReadAlignment:
    def __init__(self, name, read_seq, oriented, strand, ref_start, ref_end,
                 read_start, read_end, aligned_bp, diffs, plasmid_len,
                 blocks=None, ref_span=None):
        self.name, self.read_seq, self.oriented = name, read_seq, oriented
        self.strand = strand
        self.ref_start, self.ref_end = int(ref_start), int(ref_end)
        self.read_start, self.read_end = int(read_start), int(read_end)
        self.aligned_bp, self.diffs = int(aligned_bp), diffs
        self.plasmid_len = plasmid_len
        self.blocks = blocks
        self.ref_span = (ref_span if ref_span is not None else
                         (self.ref_end - self.ref_start) % plasmid_len + 1)
        self.quality_source = getattr(read_seq, "quality_source", "positional heuristic")

    @property
    def length(self):
        return len(self.oriented)

    @property
    def substitutions(self):
        return [d for d in self.diffs if d.kind == "sub"]

    @property
    def indels(self):
        return [d for d in self.diffs if d.kind != "sub"]

    @property
    def identity(self):
        columns = self.aligned_bp + sum(d.size for d in self.indels)
        return 100.0 * (self.aligned_bp - len(self.substitutions)) / columns if columns else 0.0

    @property
    def substitution_identity(self):
        return (100.0 * (self.aligned_bp - len(self.substitutions)) / self.aligned_bp
                if self.aligned_bp else 0.0)

    def confident_window(self, trim_start=40, trim_end=120):
        """Bounds use map-oriented read coordinates, including for reverse reads."""
        bounds = getattr(self.read_seq, "confident_bounds", None)
        lo, hi = bounds if bounds is not None else (trim_start + 1, self.length - trim_end)
        if self.strand == "-":
            lo, hi = self.length - hi + 1, self.length - lo + 1
        return lo, hi

    def ref_to_read(self, ref_pos):
        """Exact base coordinate; a deleted or unaligned base has none."""
        if self.blocks is not None:
            for r0, r1, q0, _ in self.blocks:
                offset = (ref_pos - (r0 % self.plasmid_len + 1)) % self.plasmid_len
                if offset < r1 - r0:
                    return q0 + offset + 1
            return None
        offset = (ref_pos - self.ref_start) % self.plasmid_len
        return self.read_start + offset if offset < self.ref_span else None

    def covered_positions(self):
        if self.blocks is not None:
            return {(p % self.plasmid_len) + 1 for r0, r1, _, _ in self.blocks
                    for p in range(r0, r1)}
        positions = {((self.ref_start - 1 + k) % self.plasmid_len) + 1
                     for k in range(self.ref_span)}
        for diff in self.diffs:
            if diff.kind == "del":
                positions.difference_update(diff.positions(self.plasmid_len))
        return positions

    def covers(self, start, end):
        offset = (start - self.ref_start) % self.plasmid_len
        length = (end - start) % self.plasmid_len + 1
        return offset + length <= self.ref_span


def _aligner(match=2, mismatch=-4, gap_open=-12, gap_extend=-1):
    al = Align.PairwiseAligner(mode="local")
    al.match_score, al.mismatch_score = match, mismatch
    al.open_gap_score, al.extend_gap_score = gap_open, gap_extend
    return al


def align_read(name, read_seq, pmap, aligner=None, flank=None):
    """Try both orientations against the doubled circular reference."""
    al = aligner or _aligner()
    n, target = len(pmap), pmap.doubled
    best, best_strand, best_seq = None, None, None
    for strand, sequence in (("+", str(read_seq)), ("-", revcomp(read_seq))):
        try:
            candidate = al.align(target, sequence)[0]
        except (ValueError, IndexError):
            continue
        if best is None or candidate.score > best.score:
            best, best_strand, best_seq = candidate, strand, sequence
    if best is None:
        raise ValueError(f"{name}: no alignment found")

    rb, qb = best.aligned
    blocks = [(int(r0), int(r1), int(q0), int(q1))
              for (r0, r1), (q0, q1) in zip(rb, qb)]
    diffs = []
    for r0, r1, q0, _ in blocks:
        for i in range(r1 - r0):
            x, y = best.target[r0 + i], best.query[q0 + i]
            if x != y or y not in "ACGT":
                diffs.append(Diff((r0 + i) % n + 1, q0 + i + 1, "sub", f"{x}>{y}"))
    for left, right in zip(blocks, blocks[1:]):
        dr, dq = right[0] - left[1], right[2] - left[3]
        ref_pos, read_pos = left[1] % n + 1, left[3] + 1
        if dr:
            diffs.append(Diff(ref_pos, read_pos, "del", f"-{dr} bp", dr))
        if dq:
            diffs.append(Diff(ref_pos, read_pos, "ins", f"+{dq} bp", dq))
    diffs.sort(key=lambda d: d.read)
    return ReadAlignment(
        name, read_seq, best_seq, best_strand,
        rb[0][0] % n + 1, (rb[-1][1] - 1) % n + 1,
        qb[0][0] + 1, qb[-1][1], sum(r1 - r0 for r0, r1 in rb), diffs, n,
        blocks=blocks, ref_span=int(rb[-1][1] - rb[0][0]))


def align_reads(reads, pmap, progress=True):
    al = _aligner()
    items = sorted(reads.items())
    if progress:
        try:
            from tqdm import tqdm
            items = tqdm(items, desc="aligning", unit="read")
        except ImportError:
            pass
    return [align_read(name, sequence, pmap, al) for name, sequence in items]
