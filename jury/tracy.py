"""Optional adapter for an independently installed Tracy binary.

JSON schema: https://github.com/gear-genomics/tracy/blob/main/src/json.h
The upstream basecall JSON has no trim boundaries; its quality scores are
trace-derived estimates, not a promise of calibrated Phred probabilities.
"""
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def quality_window(seq, qualities, minimum=20):
    """Longest contiguous ACGT stretch whose estimated scores meet a threshold."""
    best = (1, 0)
    start = 1
    for pos, (base, quality) in enumerate(zip(seq, qualities), 1):
        if base not in "ACGT" or quality < minimum:
            start = pos + 1
        elif pos - start > best[1] - best[0]:
            best = (start, pos)
    return best


def parse_basecalls(data):
    """Read upstream primarySeq/basecallQual without guessing other schemas."""
    if not isinstance(data, dict):
        raise ValueError("Tracy JSON must be an object")
    seq, qualities = data.get("primarySeq"), data.get("basecallQual")
    if not isinstance(seq, str) or not seq:
        raise ValueError("Tracy JSON has no nonempty primarySeq")
    seq = seq.upper()
    if set(seq) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("Tracy primarySeq contains invalid DNA characters")
    if not isinstance(qualities, list) or len(qualities) != len(seq):
        raise ValueError("Tracy basecallQual length does not match primarySeq")
    if any(type(q) is not int or not 0 <= q <= 255 for q in qualities):
        raise ValueError("Tracy basecallQual must contain integer quality scores")
    qualities = tuple(qualities)
    return seq, qualities, quality_window(seq, qualities)


def basecall(path, executable=None):
    executable = executable or shutil.which("tracy")
    if not executable:
        raise ValueError("AB1 input requires the optional Tracy executable on PATH; "
                         "install Tracy or provide FASTA/.seq basecalls")
    # Trace data is private; temporary JSON must never land in the repository.
    with tempfile.TemporaryDirectory(prefix="jury-tracy-") as temp_dir:
        output = Path(temp_dir) / "basecalls.json"
        try:
            subprocess.run([str(executable), "basecall", "-f", "json", "-o",
                            str(output), str(Path(path).resolve())],
                           check=True, capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired as exc:
            raise ValueError("Tracy basecalling timed out after 120 seconds") from exc
        except subprocess.CalledProcessError as exc:
            raise ValueError(f"Tracy basecalling failed (exit {exc.returncode})") from exc
        except OSError as exc:
            raise ValueError("could not run the Tracy executable") from exc
        try:
            data = json.loads(output.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError("Tracy did not produce valid basecall JSON") from exc
        return parse_basecalls(data)
