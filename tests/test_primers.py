import random

import pytest

from jury import primers
from jury.seqio import Map


def random_map(length=900):
    rng = random.Random(843)
    return Map("".join(rng.choice("ACGT") for _ in range(length)), [])


def covered(result):
    return {i for i, depth in enumerate(result["depth"]) if depth}


def test_reverse_primer_window_extends_toward_lower_coordinates():
    pmap = random_map(100)
    forward = primers.Primer("forward", 10, 20, "ACGT")
    reverse = primers.Primer("reverse", 80, 90, "ACGT", -1)
    result = primers.audit([forward, reverse], pmap, dead=5, good=10)
    assert covered(result) == set(range(25, 31)) | set(range(70, 76))
    assert result["windows"] == [(forward, 25, 30), (reverse, 75, 70)]


def test_reverse_primer_window_wraps_origin():
    pmap = random_map(100)
    reverse = primers.Primer("reverse", 3, 12, "ACGT", -1)
    result = primers.audit([reverse], pmap, dead=5, good=10)
    assert covered(result) == set(range(93, 99))
    assert result["pct"] == 6


def test_occurrences_counts_overlaps_and_origin_hits_once():
    pmap = Map("AAACG", [])
    assert primers.occurrences("AA", pmap) == (2, 0)
    assert primers.occurrences("GAA", pmap) == (1, 0)
    assert primers.occurrences("TTC", pmap) == (0, 1)
    assert primers.occurrences("AAACGA", pmap) == (0, 0)


def test_design_at_places_three_prime_end_inside_search_window():
    pmap = random_map()
    primer = primers.design_at(400, pmap, window=10)
    assert primer is not None
    assert 390 <= primer.three_prime < 410
    assert pmap.slice(primer.start, primer.end) == primer.seq
    assert primers.score_candidate(primer.seq, pmap) is not None


def test_fill_gaps_uses_custom_read_window_and_improves_every_step():
    pmap = random_map()
    existing = [primers.Primer("existing", 1, 20, pmap.slice(1, 20))]
    new = primers.fill_gaps(existing, pmap, spacing=200, dead=10, good=140,
                            progress=False)
    assert new
    walk = existing[:]
    missing = primers.audit(walk, pmap, dead=10, good=140)["uncovered"]
    for primer in new:
        assert primers.score_candidate(primer.seq, pmap) is not None
        walk.append(primer)
        updated = primers.audit(walk, pmap, dead=10, good=140)["uncovered"]
        assert updated < missing
        missing = updated
    assert missing == 0
    assert len({(p.start, p.end) for p in new}) == len(new)


def test_fill_gaps_reports_no_design_when_constraints_cannot_be_met():
    pmap = Map("A" * 300, [])
    new = primers.fill_gaps([], pmap, spacing=200, dead=10, good=100,
                            progress=False)
    assert new == []
    assert primers.audit(new, pmap, dead=10, good=100)["uncovered"] == 300


def test_fill_gaps_does_not_add_redundant_primers():
    pmap = random_map(100)
    existing = [primers.Primer("existing", 1, 20, pmap.slice(1, 20))]
    assert primers.fill_gaps(existing, pmap, dead=0, good=100,
                             progress=False) == []


@pytest.mark.parametrize("dead,good", [(-1, 50), (60, 50)])
def test_invalid_read_window_is_rejected(dead, good):
    with pytest.raises(ValueError, match="window"):
        primers.audit([], random_map(), dead=dead, good=good)


@pytest.mark.parametrize("spacing", [0, -1])
def test_invalid_spacing_is_rejected(spacing):
    with pytest.raises(ValueError, match="spacing"):
        primers.fill_gaps([], random_map(), spacing=spacing, progress=False)
