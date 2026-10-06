import csv
from output_utils import shorten_protein, save_orfs_to_csv


def test_shorten_protein_leaves_short_ones_alone():
    assert shorten_protein("MAKF") == "MAKF"


def test_shorten_protein_cuts_long_ones():
    assert shorten_protein("M" * 40) == "M" * 30 + "..."


def test_save_orfs_to_csv_writes_expected_columns(tmp_path):
    orfs = [{
        "sequence_id": "seq1", "orf_id": "ORF_1", "strand": "+", "frame": 1,
        "start": 1, "end": 15, "length": 15, "protein_length": 4,
        "protein": "MAKF", "dna": "ATGGCCAAATTTTAA",   # 'dna' must NOT be saved
    }]
    csv_path = tmp_path / "subfolder" / "result.csv"    # folder doesn't exist yet

    save_orfs_to_csv(orfs, csv_path)

    with open(csv_path, newline="") as csv_file:
        rows = list(csv.DictReader(csv_file))

    assert len(rows) == 1
    assert rows[0]["protein"] == "MAKF"
    assert rows[0]["start"] == "1"       # CSV stores everything as text
    assert "dna" not in rows[0]