from pathlib import Path

from sequence_utils import prepare_dna, read_fasta
from orf_finder import find_all_orfs
from translator import (
    add_proteins,
    assign_orf_ids,
    filter_by_min_length,
    find_longest_orf,
    calculate_statistics,
)
from output_utils import print_orf_table, save_orfs_to_csv

# Settings. 9 nt is tiny, chosen so our small demo sequences give results.
# For real gene finding, a typical minimum is about 300 nt (100 amino acids).
MIN_ORF_LENGTH = 9

PROJECT_FOLDER = Path(__file__).resolve().parent.parent
RESULTS_FILE = PROJECT_FOLDER / "results" / "orf_results.csv"


def get_input_sequences():
    """Ask the user for input. Returns a list of (sequence_id, clean_dna)."""
    print("Choose input type:")
    print("  1 - Type or paste a DNA sequence")
    print("  2 - Read a FASTA file")
    choice = input("Enter 1 or 2: ").strip()

    if choice == "1":
        raw_dna = input("Paste your DNA sequence: ")
        return [("manual_input", prepare_dna(raw_dna))]

    if choice == "2":
        path_text = input("FASTA file path (press Enter for data/sample.fasta): ").strip()
        if path_text == "":
            fasta_path = PROJECT_FOLDER / "data" / "sample.fasta"
        else:
            # Windows "Copy as path" adds quotation marks, so remove them
            fasta_path = Path(path_text.strip('"'))
        return read_fasta(fasta_path)

    raise ValueError("Please enter 1 or 2.")

def analyze_sequence(sequence_id, dna, min_length=MIN_ORF_LENGTH):
    """Run the ORF steps on ONE clean DNA sequence. Returns a list of ORFs."""
    orfs = find_all_orfs(dna)                         # Steps 4-5
    orfs = filter_by_min_length(orfs, MIN_ORF_LENGTH) # Step 6
    orfs = add_proteins(orfs)                         # Step 6
    orfs = assign_orf_ids(orfs)                       # Step 6

    for orf in orfs:
        orf["sequence_id"] = sequence_id              # remember which sequence it came from
    return orfs


def main():
    print("=== ORF Analysis Pipeline ===")
    print()

    try:
        sequences = get_input_sequences()
    except (ValueError, FileNotFoundError) as error:
        print("Error:", error)
        return

    all_orfs = []

    for sequence_id, dna in sequences:
        print()
        print(f"Sequence: {sequence_id} ({len(dna)} bases)")
        orfs = analyze_sequence(sequence_id, dna)

        if len(orfs) == 0:
            print(f"No ORFs of at least {MIN_ORF_LENGTH} nt found.")
            continue

        print_orf_table(orfs)

        longest = find_longest_orf(orfs)
        print(f"Longest ORF: {longest['orf_id']} ({longest['length']} nt)")

        stats = calculate_statistics(orfs)
        print(
            f"Statistics: {stats['total_orfs']} ORFs "
            f"({stats['forward_orfs']} forward, {stats['reverse_orfs']} reverse), "
            f"average length {stats['average_length']} nt"
        )

        all_orfs.extend(orfs)    # add this sequence's ORFs to the combined list

    if len(all_orfs) > 0:
        try:
            save_orfs_to_csv(all_orfs, RESULTS_FILE)
            print()
            print(f"Results saved to: {RESULTS_FILE}")
        except PermissionError:
            print()
            print("Could not save the CSV. Is it open in Excel? Close it and run again.")
    else:
        print()
        print("No ORFs found in any sequence, so no CSV was saved.")


if __name__ == "__main__":
    main()