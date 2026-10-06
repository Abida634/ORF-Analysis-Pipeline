import pytest
from sequence_utils import clean_sequence, validate_sequence, prepare_dna, read_fasta


def test_clean_removes_whitespace_and_uppercases():
    assert clean_sequence("  atg gcc aaa ") == "ATGGCCAAA"


def test_clean_handles_tabs_and_newlines():
    assert clean_sequence("AT\nGC\tAA") == "ATGCAA"


def test_validate_accepts_valid_dna():
    validate_sequence("ATGCATGC")      # should NOT raise; no assert needed


def test_validate_rejects_invalid_characters():
    with pytest.raises(ValueError, match="Invalid characters"):
        validate_sequence("ATGXC1")


def test_validate_rejects_empty():
    with pytest.raises(ValueError, match="empty"):
        validate_sequence("")


def test_prepare_dna_cleans_and_validates():
    assert prepare_dna(" atg\ngcc ") == "ATGGCC"


def test_read_fasta_joins_lines_and_cleans(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1 test record\natgc\nGGA\n")

    result = read_fasta(fasta_file)

    assert result == [("seq1", "ATGCGGA")]


def test_read_fasta_names_the_bad_record(tmp_path):
    fasta_file = tmp_path / "bad.fasta"
    fasta_file.write_text(">s1\nATGC\n>s2\nATGX\n")

    with pytest.raises(ValueError, match="Record 's2'"):
        read_fasta(fasta_file)