"""AGC-backed tests: real pulls from the HPRC archive (skipped without it)."""
from add_flanks import agc_contig, pull_regions
from project_variants import cons_agc_contig, pull_agc


def test_pull_agc_is_zero_based(agc_archive):
    # GRCh38 chr1 is N through 1-based 10000; the first base, 1-based 10001,
    # is agc position 10000. agc ranges are 0-based with an inclusive end.
    assert pull_agc(agc_archive, "GRCh38#0#chr1", 10000, 10005) == "TAACCC"


def test_pull_agc_chm13_contig(agc_archive):
    contig = cons_agc_contig("chr1")
    assert pull_agc(agc_archive, contig, 10000, 10009) == "CATATTTTAT"


def test_pull_agc_missing_contig_is_none(agc_archive):
    # agc exits 0 with no sequence for an unknown contig, so pull_agc's retry
    # loop treats it as transient; retries=1 keeps the test fast.
    assert pull_agc(agc_archive, "NOPE#0#chr1", 0, 9, retries=1) is None


def test_pull_regions_keys_by_contig(agc_archive):
    # Sample contigs (PanSN) and HG002 T2T contigs (mapped by agc_contig) come
    # back keyed by the bare contig name that add_flanks looks them up by.
    hg002 = agc_contig("HG002_1", "chr12_MATERNAL")
    regions = ["HG00097#1#CM094060.1:1000-1009", f"{hg002}:100000-100009"]
    pulled = pull_regions(agc_archive, regions, threads=1)
    assert pulled == {
        "HG00097#1#CM094060.1": "TAACCCTAAC",
        "HG002#2#chr12": "TGTTCCCTAA",
    }
