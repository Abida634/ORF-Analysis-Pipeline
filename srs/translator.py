from Bio.Seq import Seq
from orf_finder import find_all_orfs


def translate_orf(orf_dna):
    """Translate an ORF's DNA (string) into a protein (string)."""
    return str(Seq(orf_dna).translate(to_stop=True))


def add_proteins(orfs):
    """Add 'protein' and 'protein_length' to every ORF dictionary."""
    for orf in orfs:
        protein = translate_orf(orf["dna"])
        orf["protein"] = protein
        orf["protein_length"] = len(protein)
    return orfs


def filter_by_min_length(orfs, min_length):
    """Keep only ORFs whose length (in nucleotides) is at least min_length."""
    kept_orfs = []
    for orf in orfs:
        if orf["length"] >= min_length:
            kept_orfs.append(orf)
    return kept_orfs


def assign_orf_ids(orfs):
    """Give each ORF a name: ORF_1, ORF_2, ..."""
    for number, orf in enumerate(orfs, start=1):
        orf["orf_id"] = f"ORF_{number}"
    return orfs


def find_longest_orf(orfs):
    """Return the longest ORF dictionary, or None if the list is empty."""
    if len(orfs) == 0:
        return None

    longest_orf = orfs[0]
    for orf in orfs:
        if orf["length"] > longest_orf["length"]:
            longest_orf = orf
    return longest_orf


def calculate_statistics(orfs):
    """Return a dictionary of basic statistics about the ORFs."""
    if len(orfs) == 0:
        return {
            "total_orfs": 0,
            "forward_orfs": 0,
            "reverse_orfs": 0,
            "longest_length": 0,
            "shortest_length": 0,
            "average_length": 0,
        }

    lengths = [orf["length"] for orf in orfs]

    forward_count = 0
    for orf in orfs:
        if orf["strand"] == "+":
            forward_count += 1

    return {
        "total_orfs": len(orfs),
        "forward_orfs": forward_count,
        "reverse_orfs": len(orfs) - forward_count,
        "longest_length": max(lengths),
        "shortest_length": min(lengths),
        "average_length": round(sum(lengths) / len(lengths), 2),
    }


def print_orfs(orfs):
    """Print one line per ORF."""
    for orf in orfs:
        print(
            f"{orf['orf_id']} | {orf['strand']}{orf['frame']} | "
            f"start {orf['start']} | end {orf['end']} | "
            f"length {orf['length']} | protein {orf['protein']}"
        )


if __name__ == "__main__":
    # Contains: a forward ORF, a second forward ORF, and a reverse-strand ORF
    test_dna = "ATGGCCAAATTTTAACCATGAAATAGGGCTATTTCAT"

    orfs = find_all_orfs(test_dna)
    orfs = add_proteins(orfs)
    orfs = assign_orf_ids(orfs)

    print("All ORFs:")
    print_orfs(orfs)

    longest = find_longest_orf(orfs)
    print()
    print(f"Longest ORF: {longest['orf_id']} ({longest['length']} nt), "
          f"protein {longest['protein']}")

    print()
    print("Statistics:")
    for name, value in calculate_statistics(orfs).items():
        print(f"  {name}: {value}")

    print()
    print("After filtering (minimum length 12 nt):")
    filtered = filter_by_min_length(orfs, 12)
    print_orfs(filtered)