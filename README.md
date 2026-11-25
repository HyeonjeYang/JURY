# JURY

**Verify Sanger reads against a plasmid map - VERIFIED, MISMATCH, or the jury is still out.**

Check cloning results instead of Ctrl+F-ing through the map by hand.
JURY checks features, inserts, junctions and forward-strand ORFs, and plans sequencing primers.

## Install

Requires Python 3.9+. Create a virtual environment with `python -m venv .venv`.
Activate it with `.\.venv\Scripts\Activate.ps1` in PowerShell, or
`source .venv/bin/activate` on macOS/Linux, then:

```sh
python -m pip install -e .
```

Optional progress bars: `python -m pip install -e ".[progress]"`.
If `jury` is unavailable on PATH, use `python -m jury`.

## Run

Keep experimental files **outside this repository**; JURY reads them without modifying them.

- **Map / optional insert:** SnapGene `.dna`, GenBank `.gb`/`.gbk`/`.genbank`, or single-record FASTA. FASTA maps have no feature annotations.
- **Reads:** FASTA, vendor `.seq` (FASTA or plain DNA), or a directory for one clone. Matching stems and `_revComp` exports count as one reaction; Tracy-enabled AB1 takes priority.

```sh
jury verify /path/to/map.dna /path/to/reads --insert /path/to/fragment.dna --orf-motif ATGGGATCT
jury primers /path/to/map.dna --spacing 600
```

Quote paths containing spaces. For shorter vendor labels, add
`--name-regex '(HJY_\d+_\w+)-HJY'`. See `jury verify --help` and `jury primers --help`.

Maps are circular; coordinates are 1-based inclusive. Reads align in both orientations.
Identity includes substitutions and indels; deleted bases do not count as observed coverage.

| Verdict | Meaning |
|---|---|
| `VERIFIED` | Fully spanned; all spanning reads place the region inside their confidence windows, with no differences in any overlapping read. |
| `LOW-CONFIDENCE` | Fully spanned without differences, but overlaps a confidence boundary. The jury is still out. |
| `MISMATCH` | A difference or ambiguous call was observed; inspect its quality before interpreting it. |
| `PARTIAL` / `UNCOVERED` | No single read spans the region / no observed bases cover it. |

FASTA/SEQ confidence excludes the first 40 and last 120 original read bases
(`--trim-start`, `--trim-end`). This is a positional heuristic, not measured quality.
`verify` exits **1** for a feature/insert/junction/ORF difference, missing insert, or
failed requested ORF search; **2** for invalid input/options; otherwise **0**.
Zero can still mean partial, uncovered, or low-confidence: check the verdicts.

Primer planning uses annotations matching `seq_primer` (`--pattern` to change),
with reliable offsets of 50-700 bases along primer extension (`--dead`, `--good`).
Designed primers and modeled coverage still need experimental validation.

## Optional AB1 traces

Install [Tracy](https://www.gear-genomics.com/docs/tracy/installation/) separately
on the same PATH as JURY. On Windows, run both in WSL for `.ab1`/`.abi` support.
FASTA/SEQ works without Tracy; no Tracy binary or upstream source is bundled.

Tracy confidence uses the longest contiguous ACGT stretch with estimated scores >=20,
replacing positional trimming. Reports identify this source; scores are not asserted
to be calibrated Phred probabilities. Temporary JSON stays in system temp and is removed.

## Development and license

```sh
python -m pip install -e ".[dev,progress]"
python -m pytest -q
```

Tests use synthetic sequences. `.gitignore` excludes lab inputs and reports; never force-add private data.

Copyright (c) 2026 HyeonjeYang. JURY is [MIT licensed](LICENSE).
[Biopython](https://github.com/biopython/biopython/blob/master/LICENSE.rst), optional
[tqdm](https://github.com/tqdm/tqdm/blob/master/LICENCE), and optional
[Tracy](https://github.com/gear-genomics/tracy/blob/main/LICENSE) retain their own licenses.
Tracy is BSD 3-Clause, Copyright (c) 2018- European Molecular Biology Laboratory (EMBL);
no endorsement is implied. Bundled third-party code/binaries require their applicable notices.
For publications using Tracy, cite [Rausch et al., BMC Genomics 21, 230 (2020)](https://doi.org/10.1186/s12864-020-6635-8).
