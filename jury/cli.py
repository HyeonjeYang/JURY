"""Command line interface."""
import argparse
import re
import sys
from . import seqio, align, verify, primers as pr

MARK = {verify.VERIFIED: "OK ", verify.LOW_CONF: "~  ", verify.MISMATCH: "!! ",
        verify.PARTIAL: "-  ", verify.UNCOVERED: "   "}


def _nonnegative(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return number


def _positive(value):
    number = _nonnegative(value)
    if number == 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def _verify(args):
    pmap = seqio.load_map(args.map)
    reads = seqio.load_reads(args.reads, name_regex=args.name_regex)
    print(f"map   {args.map}  {len(pmap)} bp, {len(pmap.features)} features")
    print(f"reads {len(reads)}\n")

    alns = align.align_reads(reads, pmap, progress=not args.quiet)
    failed = False

    print("== reads ==")
    print(f"{'read':<18}{'len':>6}{'map span':>14}{'read span':>13}"
          f"{'ident':>8}{'sub':>5}{'indel':>6}")
    for a in alns:
        print(f"{a.name:<18}{a.length:>6}{f'{a.ref_start}-{a.ref_end}':>14}"
              f"{f'{a.read_start}-{a.read_end}':>13}{a.identity:>7.2f}%"
              f"{len(a.substitutions):>5}{len(a.indels):>6}")
        lo, hi = a.confident_window(args.trim_start, args.trim_end)
        print(f"  strand {a.strand}; confidence: {a.quality_source}; "
              f"oriented read {lo}-{hi}")

    res, depth = verify.check_features(pmap, alns, trim_start=args.trim_start,
                                       trim_end=args.trim_end, include=args.feature)
    print("\n== features ==   (confidence bounds shown per read above)")
    print(f"   {'feature':<30}{'coords':>14}{'bp':>7}{'reads':>7}{'cov':>8}"
          f"{'diffs':>7}  verdict")
    for r in res:
        print(f"{MARK[r.verdict]}{r.name[:29]:<30}{f'{r.start}-{r.end}':>14}"
              f"{r.length:>7}{r.depth_full:>7}{r.pct_covered:>7.0f}%"
              f"{len(r.diffs):>7}  {r.verdict}")

    s = verify.coverage_summary(depth, len(pmap), len(alns))
    print(f"\n== coverage ==\n  >=1 read {s['covered']}/{len(pmap)} bp "
          f"({s['pct']:.1f}%)   all {len(alns)} reads {s['pct_all']:.1f}%")
    for a, b in s["gaps"]:
        print(f"  gap {a}-{b}  ({b - a + 1} bp)")

    if args.insert:
        ins = seqio.load_map(args.insert)
        h = verify.check_insert(pmap, ins.seq)
        print(f"\n== insert ==  {args.insert}  {h['length']} bp")
        if h["found"]:
            print(f"  exact match in reference at {h['start']}-{h['end']}")
            insert_result = verify.check_region(
                "insert", h["start"], h["end"], alns, len(pmap), depth,
                args.trim_start, args.trim_end)
            print(f"  read evidence: {insert_result.verdict}; "
                  f"fully covered by {insert_result.depth_full} read(s), "
                  f"{len(insert_result.diffs)} diff(s)")
            failed |= insert_result.verdict == verify.MISMATCH
            for lab, j in (("5'", h["junction5"]), ("3'", h["junction3"])):
                junction = verify.check_region(
                    lab + " junction", j[0], j[1], alns, len(pmap), depth,
                    args.trim_start, args.trim_end)
                print(f"  {lab} junction {j[0]}/{j[1]}  "
                      f"covered by {junction.depth_full} read(s), "
                      f"{len(junction.diffs)} diff(s)  {junction.verdict}")
                failed |= junction.verdict == verify.MISMATCH
        else:
            failed = True
            print("  NOT FOUND" + (" (reverse complement is present)"
                                   if h["reverse_only"] else ""))

    if args.orf or args.orf_motif:
        try:
            orf = verify.find_orf(pmap, args.orf_motif)
            print(f"\n== ORF ==  {orf['start']}-{orf['end']}  "
                  f"{len(orf['protein'])} aa  stop {orf['stop_codon']}")
            r = verify.check_region("ORF", orf["start"], orf["end"], alns,
                                    len(pmap), depth, args.trim_start, args.trim_end)
            failed |= r.verdict == verify.MISMATCH
            print(f"  {r.verdict}  reads covering fully: {r.depth_full}, "
                  f"diffs inside: {len(r.diffs)}")
            for i in range(0, len(orf["protein"]), 60):
                print("  " + orf["protein"][i:i + 60])
        except ValueError as e:
            print(f"\n== ORF ==  {e}")
            failed = True

    bad = [r for r in res if r.verdict == verify.MISMATCH]
    by_name = {a.name: a for a in alns}
    for r in bad:
        print(f"\n!! {r.name}")
        for name, d in r.diffs:
            lo, hi = by_name[name].confident_window(args.trim_start, args.trim_end)
            confidence = ("within confidence window" if lo <= d.read <= hi
                          else "outside confidence window")
            print(f"   {name}  map {d.ref}  read {d.read}  {d.kind} {d.detail}"
                  f"  [{confidence}]")
    return 1 if bad or failed else 0


def _primers(args):
    pmap = seqio.load_map(args.map)
    existing = pr.collect_primers(pmap, args.pattern)
    print(f"map {args.map}  {len(pmap)} bp")
    print(f"\n== existing primers matching /{args.pattern}/ ==")
    if not existing:
        print("  none found")
    for p in existing:
        f, r = pr.occurrences(p.seq, pmap)
        notes = []
        if p.tm < 52:
            notes.append("low Tm")
        if len(p.seq) < 18:
            notes.append("short")
        expected = (1, 0) if p.strand >= 0 else (0, 1)
        if (f, r) != expected:
            notes.append("not unique")
        print(f"  {p.name:<24}{f'{p.start}-{p.end}':>13}{len(p.seq):>4}nt"
              f"  Tm {p.tm:5.1f}  GC {p.gc:3.0f}%  {', '.join(notes) or 'ok'}")

    a = pr.audit(existing, pmap, args.dead, args.good)
    print(f"\n== walk coverage ==  (window = {args.dead}-{args.good} bp along primer extension)")
    print(f"  {a['pct']:.1f}% covered, {len(a['gaps'])} gap(s), "
          f"{a['uncovered']} bp missing")
    for s, e in a["gaps"]:
        print(f"    gap {s}-{e}  ({e - s + 1} bp)")

    if a["gaps"] or args.force:
        new = pr.fill_gaps(existing, pmap, args.spacing, progress=not args.quiet,
                           dead=args.dead, good=args.good)
        print(f"\n== designed {len(new)} primer(s) (target spacing {args.spacing} bp) ==")
        for p in new:
            print(f"  {p.name:<12}{f'{p.start}-{p.end}':>13}{len(p.seq):>4}nt"
                  f"  Tm {p.tm:5.1f}  GC {p.gc:3.0f}%  {p.seq}")
        a2 = pr.audit(existing + new, pmap, args.dead, args.good)
        print(f"\n  combined: {len(existing) + len(new)} primers, "
              f"{a2['pct']:.1f}% covered, {len(a2['gaps'])} gap(s)")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="jury",
        description="Check Sanger reads against a SnapGene plasmid map.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    v = sub.add_parser("verify", help="align reads and judge each feature")
    v.add_argument("map", help="SnapGene .dna (or GenBank/FASTA) reference map")
    v.add_argument("reads", help="FASTA/.seq file, .ab1 (requires Tracy), or a directory")
    v.add_argument("--insert", help="map/FASTA of the cloned fragment to look for")
    v.add_argument("--orf", action="store_true", help="find and translate the longest forward-strand ORF")
    v.add_argument("--orf-motif", help="DNA motif the ORF starts with, e.g. ATGGGATCT")
    v.add_argument("--feature", nargs="*", help="only features whose name contains these")
    v.add_argument("--trim-start", type=_nonnegative, default=40,
                   help="leading read bases treated as unreliable (default 40)")
    v.add_argument("--trim-end", type=_nonnegative, default=120,
                   help="trailing read bases treated as unreliable (default 120)")
    v.add_argument("--name-regex", help=r"pull read names out of FASTA headers, e.g. '(HJY_\d+_\w+)'")
    v.add_argument("-q", "--quiet", action="store_true")
    v.set_defaults(func=_verify)

    p = sub.add_parser("primers", help="audit a sequencing primer walk and fill gaps")
    p.add_argument("map")
    p.add_argument("--pattern", default="seq_primer",
                   help="regex matching primer feature names (default seq_primer)")
    p.add_argument("--dead", type=_nonnegative, default=50,
                   help="bases after the 3' end that are not reliable (default 50)")
    p.add_argument("--good", type=_positive, default=700,
                   help="last reliable base from the 3' end (default 700)")
    p.add_argument("--spacing", type=_positive, default=600,
                   help="target 3'-end spacing when designing (default 600)")
    p.add_argument("--force", action="store_true", help="show the design step even when no gaps remain")
    p.add_argument("-q", "--quiet", action="store_true")
    p.set_defaults(func=_primers)

    args = ap.parse_args(argv)
    try:
        return args.func(args)
    except (OSError, ValueError, RuntimeError, re.error) as exc:
        ap.error(str(exc))


if __name__ == "__main__":
    sys.exit(main())
