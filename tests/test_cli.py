import random

import pytest

from jury.cli import main
from jury.seqio import Feature, Map


@pytest.fixture
def inputs(tmp_path):
    rng = random.Random(926)
    sequence = ''.join(rng.choice('ACGT') for _ in range(1000))
    reference = tmp_path / 'reference.fasta'
    reference.write_text('>reference\n' + sequence + '\n', encoding='utf-8')
    return sequence, reference, tmp_path


def test_insert_internal_mismatch_fails_without_feature_annotations(inputs, capsys):
    sequence, reference, folder = inputs
    replacement = 'A' if sequence[400] != 'A' else 'C'
    mutated = sequence[:400] + replacement + sequence[401:]
    reads = folder / 'reads.fasta'
    reads.write_text('>read\n' + mutated + '\n', encoding='utf-8')
    insert = folder / 'insert.fasta'
    insert.write_text('>insert\n' + sequence[300:500] + '\n', encoding='utf-8')
    assert main(['verify', str(reference), str(reads), '--insert', str(insert), '-q']) == 1
    output = capsys.readouterr().out
    assert 'read evidence: MISMATCH' in output
    assert "5' junction 300/301  covered by 1 read(s), 0 diff(s)" in output
    assert "3' junction 500/501  covered by 1 read(s), 0 diff(s)" in output


def test_clean_but_partial_insert_has_explicit_verdict(inputs, capsys):
    sequence, reference, folder = inputs
    reads = folder / 'reads.fasta'
    reads.write_text('>read\n' + sequence[400:700] + '\n', encoding='utf-8')
    insert = folder / 'insert.fasta'
    insert.write_text('>insert\n' + sequence[300:500] + '\n', encoding='utf-8')
    assert main(['verify', str(reference), str(reads), '--insert', str(insert), '-q']) == 0
    assert 'read evidence: PARTIAL' in capsys.readouterr().out


def test_missing_insert_and_requested_orf_fail(inputs, capsys):
    sequence, reference, folder = inputs
    insert = folder / 'missing.fasta'
    insert.write_text('>missing\n' + 'A' * 50 + '\n', encoding='utf-8')
    assert main(['verify', str(reference), str(reference), '--insert', str(insert), '-q']) == 1
    assert 'NOT FOUND' in capsys.readouterr().out
    assert main(['verify', str(reference), str(reference), '--orf-motif', 'A' * 50, '-q']) == 1
    assert 'not found' in capsys.readouterr().out


@pytest.mark.parametrize('flag,value', [('--trim-start', '-1'), ('--trim-end', '-2')])
def test_invalid_trim_rejected_before_opening_input(flag, value, capsys):
    with pytest.raises(SystemExit) as error:
        main(['verify', 'missing.dna', 'missing.seq', flag, value])
    assert error.value.code == 2
    assert 'must be nonnegative' in capsys.readouterr().err


def test_bad_input_has_concise_error(tmp_path, capsys):
    empty = tmp_path / 'empty.fasta'
    empty.write_text('', encoding='utf-8')
    with pytest.raises(SystemExit) as error:
        main(['verify', str(empty), str(empty), '-q'])
    assert error.value.code == 2
    output = capsys.readouterr().err
    assert 'error:' in output
    assert 'Traceback' not in output


def test_reverse_primer_unique_hit_is_not_flagged(monkeypatch, capsys):
    rng = random.Random(28)
    sequence = ''.join(rng.choice('ACGT') for _ in range(200))
    feature = Feature('seq_primer_reverse', 'primer_bind', 80, 100, -1)
    monkeypatch.setattr('jury.cli.seqio.load_map', lambda _: Map(sequence, [feature]))
    assert main(['primers', 'synthetic', '--dead', '0', '--good', '200', '-q']) == 0
    assert 'not unique' not in capsys.readouterr().out


def test_partial_feature_difference_reports_confidence(inputs, monkeypatch, capsys):
    sequence, _, folder = inputs
    reference_map = Map(sequence, [Feature('early feature', 'regulatory', 1, 30, 1)])
    monkeypatch.setattr('jury.cli.seqio.load_map', lambda _: reference_map)
    replacement = 'A' if sequence[20] != 'A' else 'C'
    mutated = sequence[:20] + replacement + sequence[21:]
    reads = folder / 'reads.fasta'
    reads.write_text('>read\n' + mutated[10:700] + '\n', encoding='utf-8')
    assert main(['verify', 'synthetic', str(reads), '-q']) == 1
    output = capsys.readouterr().out
    assert 'MISMATCH' in output
    assert 'outside confidence window' in output
