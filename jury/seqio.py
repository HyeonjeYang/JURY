"""Load circular maps and Sanger reads without discarding quality provenance."""
import re
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def _dna(seq, context):
    seq = str(seq).upper()
    if not seq:
        raise ValueError(f"{context}: empty DNA sequence")
    invalid = sorted(set(seq) - set("ACGTRYSWKMBDHVN"))
    if invalid:
        raise ValueError(f"{context}: invalid DNA characters: {''.join(invalid)!r}")
    return seq


class Feature:
    __slots__ = ("name", "type", "start", "end", "strand", "plasmid_len")

    def __init__(self, name, type_, start, end, strand, plasmid_len=None):
        self.name, self.type = name, type_
        self.start, self.end = start, end   # 1-based inclusive
        self.strand = strand
        self.plasmid_len = plasmid_len

    @property
    def length(self):
        if self.end < self.start and self.plasmid_len:
            return (self.end - self.start) % self.plasmid_len + 1
        return self.end - self.start + 1

    def __repr__(self):
        return f"<{self.type} {self.name} {self.start}-{self.end}{'+' if self.strand >= 0 else '-'}>"


class Map:
    """A plasmid map: sequence plus annotated features."""

    def __init__(self, seq, features, name=""):
        self.seq = _dna(seq, "map")
        self.features = features
        self.name = name

    def __len__(self):
        return len(self.seq)

    @property
    def doubled(self):
        """Sequence repeated once, so circular spans can be sliced linearly."""
        return self.seq + self.seq

    def find(self, motif):
        """1-based position of motif, or -1."""
        motif = _dna(motif, "motif")
        if len(motif) > len(self):
            return -1
        i = self.doubled.find(motif)
        return -1 if i < 0 else i % len(self) + 1

    def slice(self, start, end):
        """1-based inclusive slice that may wrap the origin."""
        n = len(self)
        length = end - start + 1 if end >= start else (end - start) % n + 1
        offset = (start - 1) % n
        copies = (offset + length + n - 1) // n
        return (self.seq * copies)[offset:offset + length]

    def features_named(self, *patterns):
        pats = [re.compile(p, re.I) for p in patterns]
        return [f for f in self.features if any(p.search(f.name) for p in pats)]


def _feature_bounds(location, n, name):
    intervals = sorted((int(part.start), int(part.end)) for part in location.parts)
    merged = []
    for start, end in intervals:
        if not 0 <= start < end <= n:
            raise ValueError(f"feature {name!r} has invalid coordinates")
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    if len(merged) == 1:
        return merged[0][0] + 1, merged[0][1]
    if len(merged) == 2 and merged[0][0] == 0 and merged[-1][1] == n:
        return merged[-1][0] + 1, merged[0][1]
    raise ValueError(f"feature {name!r} has disjoint segments; only contiguous circular features are supported")


def load_map(path):
    """Read a SnapGene .dna file, or a FASTA/GenBank fallback."""
    path = Path(path)
    fmt = {".dna": "snapgene", ".gb": "genbank", ".gbk": "genbank",
           ".genbank": "genbank"}.get(path.suffix.lower(), "fasta")
    rec = SeqIO.read(str(path), fmt)
    feats = []
    for f in rec.features:
        if f.type == "source":
            continue
        q = f.qualifiers
        name = (q.get("label") or q.get("name") or q.get("gene")
                or q.get("product") or ["?"])[0]
        if f.location is None:
            raise ValueError(f"feature {name!r} has no usable location")
        start, end = _feature_bounds(f.location, len(rec), name)
        feats.append(Feature(name, f.type, start, end,
                             1 if f.location.strand in (1, None) else -1,
                             len(rec)))
    return Map(rec.seq, feats, path.stem)


class Read(str):
    """String-compatible read; bounds use the original read orientation."""

    def __new__(cls, seq, qualities=None, confident_bounds=None,
                quality_source="positional heuristic"):
        obj = super().__new__(cls, _dna(seq, "read"))
        obj.qualities = tuple(qualities) if qualities is not None else None
        if obj.qualities is not None and len(obj.qualities) != len(obj):
            raise ValueError("read sequence and quality lengths differ")
        obj.confident_bounds = confident_bounds
        obj.quality_source = quality_source
        return obj


_TEXT_READS = {".fa", ".fasta", ".fas", ".seq", ".fna"}
_TRACE_READS = {".ab1", ".abi"}


def _sample_key(path):
    return re.sub(r"_revcomp$", "", path.stem, flags=re.I).casefold()


def _read_files(path, tracy_executable):
    if not path.is_dir():
        return [path]
    groups = {}
    for fp in sorted(path.iterdir()):
        if fp.is_file() and fp.suffix.lower() in _TEXT_READS | _TRACE_READS:
            groups.setdefault(_sample_key(fp), []).append(fp)
    files = []
    for group in groups.values():
        usable = [fp for fp in group if fp.suffix.lower() not in _TRACE_READS or tracy_executable]
        if not usable:
            raise ValueError("AB1 input requires the optional Tracy executable on PATH; "
                             "install Tracy or provide FASTA/.seq basecalls")

        def preference(fp):
            trace = fp.suffix.lower() in _TRACE_READS
            reverse = fp.stem.casefold().endswith("_revcomp")
            return (0 if trace else 2 if reverse else 1,
                    0 if fp.suffix.lower() == ".seq" else 1, fp.name)

        files.append(min(usable, key=preference))
    return sorted(files)


def load_reads(path, name_regex=None):
    r"""Load FASTA/.seq or Tracy-basecalled AB1 input as {name: Read}.

    Directory inputs group same-stem files (including _revComp exports) as one
    reaction. Prefer AB1 with Tracy, otherwise original text over reverse copies.
    name_regex: optional pattern; its first group becomes the read name.
    """
    path = Path(path)
    tracy_executable = shutil.which("tracy")
    files = _read_files(path, tracy_executable)
    if not files:
        raise FileNotFoundError(f"no FASTA/.seq/AB1 files under {path}")
    pat = re.compile(name_regex) if name_regex else None
    reads = {}
    for fp in files:
        if fp.suffix.lower() in _TRACE_READS:
            from .tracy import basecall
            seq, qualities, bounds = basecall(fp, executable=tracy_executable)
            records = [(fp.stem, fp.stem, Read(seq, qualities, bounds,
                        "tracy estimated trace quality (score >=20)"))]
        elif fp.suffix.lower() in _TEXT_READS:
            text = fp.read_text(encoding="utf-8-sig")
            if not text.strip():
                raise ValueError(f"{fp.name}: no reads found")
            if not text.lstrip().startswith(">"):
                if fp.suffix.lower() != ".seq":
                    raise ValueError(f"{fp.name}: expected FASTA headers starting with '>'")
                records = [(fp.stem, fp.stem, Read("".join(text.split())))]
            else:
                from io import StringIO
                records = [(rec.id, rec.description, Read(rec.seq))
                           for rec in SeqIO.parse(StringIO(text.lstrip()), "fasta")]
        else:
            raise ValueError(f"unsupported read format: {fp.suffix or '(no extension)'}")
        if not records:
            raise ValueError(f"{fp.name}: no reads found")
        for name, description, seq in records:
            if pat:
                m = pat.search(description) or pat.search(fp.name)
                if m:
                    name = m.group(1) if m.groups() else m.group(0)
            while name in reads:
                name += "'"
            reads[name] = seq
    return reads


def revcomp(s):
    return str(Seq(s).reverse_complement())
