from orf_finder import (
    get_reverse_complement,
    get_codons,
    find_orfs_in_frame,
    find_all_orfs,
)


def test_reverse_complement_basic():
    assert get_reverse_complement("ATGGCC") == "GGCCAT"


def test_reverse_complement_twice_gives_original():
    dna = "ATGGCCAAATTTTAA"
    assert get_reverse_complement(get_reverse_complement(dna)) == dna


def test_get_codons_frame_1():
    assert get_codons("ATGGCCAAATTTTAA", 1) == ["ATG", "GCC", "AAA", "TTT", "TAA"]


def test_get_codons_frame_2_ignores_leftover():
    # Frame 2 starts at the 2nd base; the last "AA" is an incomplete codon
    assert get_codons("ATGGCCAAATTTTAA", 2) == ["TGG", "CCA", "AAT", "TTT"]


def test_find_orfs_in_frame_simple():
    assert find_orfs_in_frame(["CCC", "ATG", "AAA", "TAG"]) == [(1, 3)]


def test_find_orfs_in_frame_nested_atg_uses_first():
    assert find_orfs_in_frame(["ATG", "ATG", "AAA", "TAG"]) == [(0, 3)]


def test_find_orfs_in_frame_no_stop_gives_nothing():
    assert find_orfs_in_frame(["ATG", "GCC", "AAA"]) == []


def test_find_orfs_in_frame_two_orfs():
    codons = ["ATG", "TAA", "ATG", "AAA", "TGA"]
    assert find_orfs_in_frame(codons) == [(0, 1), (2, 4)]


def test_find_all_orfs_forward():
    orfs = find_all_orfs("ATGGCCAAATTTTAA")

    assert len(orfs) == 1
    assert orfs[0]["strand"] == "+"
    assert orfs[0]["frame"] == 1
    assert orfs[0]["start"] == 1
    assert orfs[0]["end"] == 15
    assert orfs[0]["length"] == 15


def test_find_all_orfs_reverse():
    orfs = find_all_orfs("GGCTATTTCATCC")

    assert len(orfs) == 1
    assert orfs[0]["strand"] == "-"
    assert orfs[0]["frame"] == 3
    assert orfs[0]["start"] == 3
    assert orfs[0]["end"] == 11
    assert orfs[0]["dna"] == "ATGAAATAG"


def test_reverse_orf_coordinates_map_back_to_original():
    # The region at start..end of the ORIGINAL sequence must be
    # the reverse complement of the ORF's DNA.
    dna = "GGCTATTTCATCC"
    orf = find_all_orfs(dna)[0]

    region = dna[orf["start"] - 1 : orf["end"]]

    assert get_reverse_complement(region) == orf["dna"]


def test_find_all_orfs_no_orfs():
    assert find_all_orfs("CCCCCC") == []