import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation

from jury.seqio import Map, Read, load_map, load_reads


def test_circular_slice_normalizes_unfolded_coordinates():
    pmap = Map("ACGTAC", [])
    assert pmap.slice(5, 2) == "ACAC"
    assert pmap.slice(11, 14) == "ACAC"
    assert pmap.slice(-1, 2) == "ACAC"
    assert pmap.find("CA") == 6
    assert pmap.find("ACGTACA") == -1


@pytest.mark.parametrize("strand", [1, -1])
def test_origin_spanning_genbank_feature_retains_circular_interval(tmp_path, strand):
    record = SeqRecord(Seq("ACGT" * 25), id="synthetic")
    record.annotations["molecule_type"] = "DNA"
    parts = [SimpleLocation(90, 100, strand), SimpleLocation(0, 10, strand)]
    if strand < 0:
        parts.reverse()
    record.features = [SeqFeature(CompoundLocation(parts), type="CDS",
                                  qualifiers={"label": ["origin"]})]
    path = tmp_path / "map.gb"
    SeqIO.write(record, path, "genbank")
    feature = load_map(path).features[0]
    assert (feature.start, feature.end, feature.length, feature.strand) == (91, 10, 20, strand)


def test_disjoint_annotation_is_not_misreported_as_one_interval(tmp_path):
    record = SeqRecord(Seq("ACGT" * 25), id="synthetic")
    record.annotations["molecule_type"] = "DNA"
    record.features = [SeqFeature(CompoundLocation([SimpleLocation(5, 10),
                                                   SimpleLocation(20, 30)]), type="CDS")]
    path = tmp_path / "map.gb"
    SeqIO.write(record, path, "genbank")
    with pytest.raises(ValueError, match="disjoint"):
        load_map(path)


def test_vendor_exports_count_once_without_tracy(tmp_path, monkeypatch):
    monkeypatch.setattr("jury.seqio.shutil.which", lambda name: None)
    (tmp_path / "sample.seq").write_text(">vendor-sample\nACGTAC\n", encoding="utf-8")
    (tmp_path / "sample_revComp.fasta").write_text(">reverse\nGTACGT\n", encoding="utf-8")
    (tmp_path / "sample.ab1").write_bytes(b"synthetic placeholder")
    reads = load_reads(tmp_path)
    assert list(reads) == ["vendor-sample"]
    assert reads["vendor-sample"] == "ACGTAC"
    assert reads["vendor-sample"].quality_source == "positional heuristic"
    assert reads["vendor-sample"].confident_bounds is None


def test_vendor_exports_prefer_tracy_when_available(tmp_path, monkeypatch):
    monkeypatch.setattr("jury.seqio.shutil.which", lambda name: "tracy")
    monkeypatch.setattr("jury.tracy.basecall", lambda *a, **kw: ("ACGT", (10, 30, 30, 10), (2, 3)))
    (tmp_path / "sample.seq").write_text(">sample\nAAAA\n", encoding="utf-8")
    (tmp_path / "sample.ab1").write_bytes(b"synthetic placeholder")
    reads = load_reads(tmp_path)
    assert len(reads) == 1
    assert reads["sample"] == "ACGT"
    assert reads["sample"].confident_bounds == (2, 3)
    assert reads["sample"].qualities == (10, 30, 30, 10)
    assert "tracy" in reads["sample"].quality_source


def test_plain_vendor_seq_and_regex(tmp_path):
    path = tmp_path / "HJY_1_alpha.seq"
    path.write_text("acgt\nacgt\n", encoding="utf-8")
    reads = load_reads(path, r"(HJY_\d+)")
    assert list(reads) == ["HJY_1"]
    assert reads["HJY_1"] == "ACGTACGT"


@pytest.mark.parametrize("content", ["", ">empty\n", ">bad\nACGT!\n"])
def test_invalid_reads_fail_clearly(tmp_path, content):
    path = tmp_path / "bad.fasta"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError):
        load_reads(path)


def test_trace_without_tracy_has_actionable_error(tmp_path, monkeypatch):
    monkeypatch.setattr("jury.seqio.shutil.which", lambda name: None)
    path = tmp_path / "sample.ab1"
    path.write_bytes(b"synthetic placeholder")
    with pytest.raises(ValueError, match="optional Tracy"):
        load_reads(path)


def test_empty_map_rejected():
    with pytest.raises(ValueError, match="empty DNA"):
        Map("", [])


def test_string_compatible_read_preserves_quality_metadata():
    read = Read("acgt", [30, 30, 30, 30], (1, 4), "synthetic quality")
    assert isinstance(read, str)
    assert read == "ACGT"
    assert read.qualities == (30, 30, 30, 30)
