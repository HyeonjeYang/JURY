"""Turn read alignments into per-feature verdicts."""
from Bio.Seq import Seq
from .seqio import revcomp

VERIFIED = "VERIFIED"          # covered inside the confident window, no diffs
LOW_CONF = "LOW-CONFIDENCE"    # covered, but overlaps a read's ramp-up or tail
MISMATCH = "MISMATCH"          # diffs found inside the feature
PARTIAL = "PARTIAL"            # only part of it is covered
UNCOVERED = "UNCOVERED"


class RegionResult:
    __slots__ = ("name", "start", "end", "depth_min", "depth_full", "pct_covered",
                 "diffs", "in_window", "verdict", "read_spans", "plasmid_len")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    @property
    def length(self):
        return (self.end - self.start) % self.plasmid_len + 1


def _positions(start, end, n):
    return [((start - 1 + k) % n) + 1 for k in range((end - start) % n + 1)]


def depth_profile(alns, plasmid_len):
    """Count observed bases; a deletion is evidence, but not a basecall."""
    depth = [0] * (plasmid_len + 1)
    for alignment in alns:
        for position in alignment.covered_positions():
            depth[position] += 1
    return depth


def check_region(name, start, end, alns, plasmid_len, depth,
                 trim_start=40, trim_end=120):
    """Judge an interval using all overlapping evidence, including partial reads."""
    positions = set(_positions(start, end, plasmid_len))
    covered = [p for p in positions if depth[p] > 0]
    pct = 100.0 * len(covered) / len(positions)
    dmin = min((depth[p] for p in positions), default=0)
    full = [a for a in alns if a.covers(start, end)]
    diffs, spans = [], {}
    in_window = bool(full)
    for alignment in alns:
        for diff in alignment.diffs:
            if positions.intersection(diff.positions(plasmid_len)):
                diffs.append((alignment.name, diff))
    for alignment in full:
        lo, hi = alignment.confident_window(trim_start, trim_end)
        rs, re_ = alignment.ref_to_read(start), alignment.ref_to_read(end)
        spans[alignment.name] = (rs, re_)
        if rs is None or re_ is None or rs < lo or re_ > hi:
            in_window = False
    if diffs:
        verdict = MISMATCH
    elif not full:
        verdict = UNCOVERED if pct == 0 else PARTIAL
    elif not in_window:
        verdict = LOW_CONF
    else:
        verdict = VERIFIED
    return RegionResult(name=name, start=start, end=end, depth_min=dmin,
                        depth_full=len(full), pct_covered=pct, diffs=diffs,
                        in_window=in_window, verdict=verdict, read_spans=spans,
                        plasmid_len=plasmid_len)


def check_features(pmap, alns, types=("CDS", "promoter", "RBS", "terminator",
                                      "protein_bind", "misc_feature", "regulatory"),
                   trim_start=40, trim_end=120, include=None):
    """Judge every annotated feature of interest."""
    depth = depth_profile(alns, len(pmap))
    out = []
    for f in pmap.features:
        if include and not any(s.lower() in f.name.lower() for s in include):
            continue
        if not include and f.type not in types:
            continue
        out.append(check_region(f.name, f.start, f.end, alns, len(pmap), depth,
                                trim_start, trim_end))
    out.sort(key=lambda r: r.start)
    return out, depth


def coverage_summary(depth, plasmid_len, n_reads):
    """Plasmid-wide coverage and the uncovered stretches."""
    cov1 = sum(1 for p in range(1, plasmid_len + 1) if depth[p] >= 1)
    covall = sum(1 for p in range(1, plasmid_len + 1) if n_reads > 0 and depth[p] >= n_reads)
    gaps, run = [], None
    for p in range(1, plasmid_len + 1):
        if depth[p] == 0 and run is None:
            run = p
        elif depth[p] > 0 and run is not None:
            gaps.append((run, p - 1)); run = None
    if run is not None:
        gaps.append((run, plasmid_len))
    return {"covered": cov1, "pct": 100.0 * cov1 / plasmid_len,
            "all_reads": covall, "pct_all": 100.0 * covall / plasmid_len,
            "gaps": gaps, "uncovered": plasmid_len - cov1}


def find_orf(pmap, motif=None, min_aa=50):
    """Locate an ORF: from a start motif if given, else the longest one."""
    n = len(pmap)
    if motif:
        s = pmap.find(motif)
        if s < 0:
            raise ValueError(f"start motif {motif!r} not found in map")
        return _extend(pmap, s)
    best = None
    for s in range(1, n + 1):
        if pmap.slice(s, s + 2) != "ATG":
            continue
        try:
            orf = _extend(pmap, s)
        except ValueError:
            continue
        if len(orf["protein"]) >= min_aa and (best is None or
                                              len(orf["protein"]) > len(best["protein"])):
            best = orf
    if best is None:
        raise ValueError("no ORF found")
    return best


def _extend(pmap, start):
    n = len(pmap)
    codons = []
    for k in range(0, n, 3):
        c = pmap.slice(start + k, start + k + 2)
        if len(c) < 3:
            break
        if c in ("TAA", "TAG", "TGA"):
            prot = str(Seq("".join(codons)).translate())
            return {"start": start, "stop_start": (start + k - 1) % n + 1,
                    "end": (start + k + 1) % n + 1, "protein": prot,
                    "stop_codon": c, "nt": "".join(codons)}
        codons.append(c)
    raise ValueError("no in-frame stop codon")


def check_insert(pmap, insert_seq):
    """Is the insert present verbatim, and where are its junctions?"""
    ins = str(insert_seq).upper()
    n = len(pmap)
    if not ins or len(ins) > n:
        raise ValueError("insert must contain between 1 and one map length of bases")
    i = pmap.doubled.find(ins)
    rc_hit = pmap.doubled.find(revcomp(ins))
    if i < 0:
        return {"found": False, "reverse_only": rc_hit >= 0, "length": len(ins)}
    start = i % n + 1
    end = (i + len(ins) - 1) % n + 1
    return {"found": True, "reverse_only": False, "length": len(ins),
            "start": start, "end": end,
            "junction5": ((start - 2) % n + 1, start),
            "junction3": (end, end % n + 1)}
