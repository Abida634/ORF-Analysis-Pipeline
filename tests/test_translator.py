from translator import (
    translate_orf,
    add_proteins,
    filter_by_min_length,
    assign_orf_ids,
    find_longest_orf,
    calculate_statistics,
)


def test_translate_orf_simple():
    assert translate_orf("ATGAAATAG") == "MK"


def test_translate_orf_has_no_stop_symbol():
    protein = translate_orf("ATGGCCAAATTTTAA")
    assert protein == "MAKF"
    assert "*" not in protein


def test_add_proteins_adds_protein_and_length():
    orfs = [{"dna": "ATGAAATAG"}]
    result = add_proteins(orfs)

    assert result[0]["protein"] == "MK"
    assert result[0]["protein_length"] == 2


def test_filter_by_min_length_keeps_equal_and_longer():
    orfs = [{"length": 6}, {"length": 9}, {"length": 15}]
    kept = filter_by_min_length(orfs, 9)

    assert [orf["length"] for orf in kept] == [9, 15]


def test_assign_orf_ids_starts_at_one():
    orfs = [{}, {}]
    assign_orf_ids(orfs)

    assert orfs[0]["orf_id"] == "ORF_1"
    assert orfs[1]["orf_id"] == "ORF_2"


def test_find_longest_orf():
    orfs = [{"length": 9}, {"length": 15}, {"length": 12}]
    assert find_longest_orf(orfs)["length"] == 15


def test_find_longest_orf_empty_list_gives_none():
    assert find_longest_orf([]) is None


def test_calculate_statistics():
    orfs = [
        {"length": 15, "strand": "+"},
        {"length": 9, "strand": "+"},
        {"length": 9, "strand": "-"},
    ]
    stats = calculate_statistics(orfs)

    assert stats["total_orfs"] == 3
    assert stats["forward_orfs"] == 2
    assert stats["reverse_orfs"] == 1
    assert stats["longest_length"] == 15
    assert stats["shortest_length"] == 9
    assert stats["average_length"] == 11.0


def test_calculate_statistics_empty_list():
    stats = calculate_statistics([])

    assert stats["total_orfs"] == 0
    assert stats["average_length"] == 0