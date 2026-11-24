"""Regression cases generated from synthetic DNA; no laboratory data is stored here."""
import random

import pytest

from jury.align import Diff, ReadAlignment, align_read
from jury.seqio import Map, revcomp
from jury.verify import (
    LOW_CONF,
    MISMATCH,
    PARTIAL,
    UNCOVERED,
    VERIFIED,
    check_region,
    coverage_summary,
    depth_profile,
)


def random_map(length=1000):
    rng = random.Random(621)
    return Map("".join(rng.choice("ACGT") for _ in range(length)), [])


def linear_alignment(name, start, end, diffs=(), plasmid_len=1000):
    length = end - start + 1
    sequence = "A" * length
    return ReadAlignment(name, sequence, sequence, "+", start, end, 1, length,
                         length, list(diffs), plasmid_len)


def region(start, end, alignments, plasmid_len=1000, **kwargs):
    depth = depth_profile(alignments, plasmid_len)
    return check_region("synthetic region", start, end, alignments, plasmid_len,
                        depth, **kwargs)


def test_partial_read_can_contradict_an_otherwise_verified_feature():
    clean = linear_alignment("clean", 100, 900)
    partial = linear_alignment("partial", 450, 700,
                               [Diff(500, 51, "sub", "A>G")])
    result = region(400, 600, [clean, partial], trim_start=0, trim_end=0)
    assert result.verdict == MISMATCH
    assert result.depth_full == 1
    assert [(name, diff.ref) for name, diff in result.diffs] == [("partial", 500)]


def test_mutation_in_a_partial_read_is_reported_without_a_full_read():
    alignment = linear_alignment("partial", 450, 700,
                                 [Diff(500, 51, "sub", "A>G")])
    result = region(400, 600, [alignment], trim_start=0, trim_end=0)
    assert result.verdict == MISMATCH
    assert result.depth_full == 0
    assert result.pct_covered == pytest.approx(100 * 151 / 201)


def test_region_entirely_inside_a_deletion_is_not_verified_or_covered():
    # Map 100..489 and 520..900 have observed bases; 490..519 were deleted.
    sequence = "A" * 771
    alignment = ReadAlignment(
        "deleted", sequence, sequence, "+", 100, 900, 1, 771, 771,
        [Diff(490, 391, "del", "-30 bp", size=30)], 1000,
        blocks=[(99, 489, 0, 390), (519, 900, 390, 771)], ref_span=801,
    )
    result = region(500, 510, [alignment], trim_start=0, trim_end=0)
    assert result.verdict == MISMATCH
    assert result.pct_covered == 0
    assert result.depth_min == 0
    assert alignment.ref_to_read(489) == 390
    assert alignment.ref_to_read(500) is None
    assert alignment.ref_to_read(520) == 391
    assert alignment.ref_to_read(901) is None


def test_deletion_evidence_and_missing_coverage_wrap_the_origin():
    sequence = "A" * 25
    deletion = Diff(98, 9, "del", "-6 bp", size=6)
    alignment = ReadAlignment(
        "origin deletion", sequence, sequence, "+", 90, 20, 1, 25, 25,
        [deletion], 100,
        blocks=[(89, 97, 0, 8), (103, 120, 8, 25)], ref_span=31,
    )
    assert deletion.positions(100) == {98, 99, 100, 1, 2, 3}
    result = region(99, 2, [alignment], plasmid_len=100, trim_start=0, trim_end=0)
    assert result.verdict == MISMATCH
    assert result.length == 4
    assert result.pct_covered == 0
    assert alignment.ref_to_read(97) == 8
    assert alignment.ref_to_read(4) == 9
    assert alignment.ref_to_read(1) is None


def test_insertions_shift_read_coordinates_without_removing_reference_coverage():
    sequence = "A" * 306
    alignment = ReadAlignment(
        "inserted", sequence, sequence, "+", 100, 400, 1, 306, 301,
        [Diff(201, 102, "ins", "+5 bp", size=5)], 1000,
        blocks=[(99, 200, 0, 101), (200, 400, 106, 306)], ref_span=301,
    )
    assert alignment.ref_to_read(200) == 101
    assert alignment.ref_to_read(201) == 107
    assert alignment.ref_to_read(400) == 306
    assert alignment.covered_positions() == set(range(100, 401))
    result = region(190, 220, [alignment], trim_start=0, trim_end=0)
    assert result.verdict == MISMATCH
    assert result.pct_covered == 100
    assert alignment.identity == pytest.approx(100 * 301 / 306)
    assert alignment.substitution_identity == 100


def test_identity_accounts_for_deleted_bases_and_substitutions():
    alignment = linear_alignment("differences", 1, 100,
                                 [Diff(20, 20, "sub", "A>G"),
                                  Diff(60, 60, "del", "-3 bp", size=3)])
    alignment.aligned_bp = 97
    assert alignment.identity == pytest.approx(96)
    assert alignment.substitution_identity == pytest.approx(100 * 96 / 97)


@pytest.mark.parametrize("strand,expected", [("+", (41, 880)), ("-", (121, 960))])
def test_heuristic_confidence_uses_original_read_length_and_orientation(strand, expected):
    # Local alignment clips both ends; clipping must not apply an extra tail trim.
    sequence = "A" * 1000
    alignment = ReadAlignment("clipped", sequence, sequence, strand,
                              100, 800, 100, 800, 701, [], 2000)
    assert alignment.confident_window() == expected


@pytest.mark.parametrize("strand,expected", [("+", (21, 950)), ("-", (51, 980))])
def test_measured_confidence_bounds_override_heuristic_and_follow_orientation(strand, expected):
    class QualityRead(str):
        confident_bounds = (21, 950)
        quality_source = "synthetic measured quality"

    read = QualityRead("A" * 1000)
    alignment = ReadAlignment("measured", read, str(read), strand,
                              1, 1000, 1, 1000, 1000, [], 2000)
    assert alignment.confident_window(trim_start=100, trim_end=200) == expected
    assert alignment.quality_source == "synthetic measured quality"


def test_exact_bases_in_the_ramp_up_remain_low_confidence():
    alignment = linear_alignment("ramp up", 1, 1000)
    assert region(10, 30, [alignment]).verdict == LOW_CONF
    assert region(100, 200, [alignment]).verdict == VERIFIED
    assert region(900, 950, [alignment]).verdict == LOW_CONF


def test_partial_and_absent_coverage_keep_distinct_verdicts():
    alignment = linear_alignment("partial", 100, 200)
    assert region(50, 150, [alignment]).verdict == PARTIAL
    assert region(300, 350, [alignment]).verdict == UNCOVERED


def test_long_origin_crossing_read_is_not_truncated_by_a_fixed_flank():
    pmap = random_map(3200)
    read = pmap.slice(3001, 1800)
    alignment = align_read("long circular read", read, pmap)
    assert len(read) == alignment.aligned_bp == 2000
    assert (alignment.ref_start, alignment.ref_end) == (3001, 1800)
    assert alignment.identity == 100
    assert alignment.covers(3100, 100)
    assert not alignment.covers(1700, 3100)
    assert alignment.ref_to_read(3200) == 200
    assert alignment.ref_to_read(1) == 201
    assert region(3100, 100, [alignment], plasmid_len=3200).verdict == VERIFIED


def test_reverse_complement_read_is_aligned_to_map_coordinates():
    pmap = random_map()
    read = revcomp(pmap.slice(201, 800))
    alignment = align_read("reverse", read, pmap)
    assert alignment.strand == "-"
    assert (alignment.ref_start, alignment.ref_end) == (201, 800)
    assert alignment.aligned_bp == 600
    assert alignment.identity == 100
    assert alignment.confident_window() == (121, 560)


def test_alignment_retains_multibase_deletion_evidence():
    pmap = random_map()
    read = pmap.seq[99:489] + pmap.seq[519:900]
    alignment = align_read("synthetic deletion", read, pmap)
    assert [(diff.kind, diff.size) for diff in alignment.diffs] == [("del", 30)]
    assert (alignment.ref_start, alignment.ref_end) == (100, 900)
    assert alignment.ref_to_read(480) == 381
    assert alignment.ref_to_read(500) is None
    assert alignment.ref_to_read(520) == 391
    assert alignment.identity == pytest.approx(100 * 771 / 801)
    assert region(500, 510, [alignment]).verdict == MISMATCH


def test_alignment_retains_multibase_insertion_evidence():
    pmap = random_map()
    read = pmap.seq[99:200] + "ACGTA" + pmap.seq[200:400]
    alignment = align_read("synthetic insertion", read, pmap)
    assert [(diff.kind, diff.size) for diff in alignment.diffs] == [("ins", 5)]
    assert alignment.ref_to_read(200) == 101
    assert alignment.ref_to_read(201) == 107
    assert alignment.identity == pytest.approx(100 * 301 / 306)
    assert region(190, 220, [alignment]).verdict == MISMATCH


def test_matching_ambiguous_bases_do_not_count_as_verified_sequence():
    pmap = random_map()
    sequence = pmap.seq[:250] + "N" + pmap.seq[251:]
    ambiguous_map = Map(sequence, [])
    alignment = align_read("ambiguous", ambiguous_map.slice(101, 400), ambiguous_map)
    result = region(240, 260, [alignment], trim_start=0, trim_end=0)
    assert result.verdict == MISMATCH
    assert any(diff.ref == 251 and diff.kind == "sub" for _, diff in result.diffs)


def test_no_reads_never_implies_complete_all_read_coverage():
    summary = coverage_summary(depth_profile([], 100), 100, 0)
    assert summary["covered"] == summary["all_reads"] == 0
    assert summary["pct"] == summary["pct_all"] == 0
    assert summary["gaps"] == [(1, 100)]
