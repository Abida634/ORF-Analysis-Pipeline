import csv
from pathlib import Path

# The columns of the CSV file, in this order
CSV_COLUMNS = [
    "sequence_id", "orf_id", "strand", "frame",
    "start", "end", "length", "protein_length", "protein",
]


def shorten_protein(protein, max_characters=30):
    """Shorten long proteins for terminal display (the CSV keeps the full one)."""
    if len(protein) > max_characters:
        return protein[:max_characters] + "..."
    return protein


def print_orf_table(orfs):
    """Print the ORFs as an aligned table."""
    header = (
        f"{'ORF':<8} {'Strand':<7} {'Frame':<6} "
        f"{'Start':>6} {'End':>6} {'Length':>7}  Protein"
    )
    print(header)
    print("-" * len(header))

    for orf in orfs:
        print(
            f"{orf['orf_id']:<8} {orf['strand']:<7} {orf['frame']:<6} "
            f"{orf['start']:>6} {orf['end']:>6} {orf['length']:>7}  "
            f"{shorten_protein(orf['protein'])}"
        )


def save_orfs_to_csv(orfs, csv_path):
    """Save the ORF dictionaries to a CSV file."""
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)   # make sure results/ exists

    with open(csv_path, "w", newline="") as csv_file:
        writer = csv.DictWriter(
            csv_file, fieldnames=CSV_COLUMNS, extrasaction="ignore"
        )
        writer.writeheader()      # first line: column names
        writer.writerows(orfs)    # one line per ORF dictionary