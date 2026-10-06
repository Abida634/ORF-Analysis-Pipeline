from Bio.Seq import Seq
START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}

def get_reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def get_codons(sequence, frame_number):

    offset = frame_number - 1
    codons = []

    for position in range(offset, len(sequence) - 2, 3):
        codon = sequence[position:position + 3]
        codons.append(codon)

    return codons


def get_six_frames(sequence):

    reverse_sequence = get_reverse_complement(sequence)

    strands = [("+", sequence), ("-", reverse_sequence)]
    all_frames = []

    for strand_symbol, strand_sequence in strands:
        for frame_number in [1, 2, 3]:
            codons = get_codons(strand_sequence, frame_number)
            all_frames.append((strand_symbol, frame_number, codons))

    return all_frames

def find_orfs_in_frame(codons):
    """Scan one frame's codon list.

    Returns a list of (start_codon_index, stop_codon_index) pairs,
    one pair for each complete ORF (ATG ... stop).
    """
    orf_positions = []
    start_index = None            # None means "not inside an ORF right now"

    for index, codon in enumerate(codons):
        if start_index is None:
            # Looking for a start codon
            if codon == START_CODON:
                start_index = index
        else:
            # Inside an ORF, looking for a stop codon
            if codon in STOP_CODONS:
                orf_positions.append((start_index, index))
                start_index = None    # ORF closed; go back to searching

    return orf_positions


def find_all_orfs(sequence):
    """Find ORFs in all 6 frames. Returns a list of dictionaries."""
    sequence_length = len(sequence)
    reverse_sequence = get_reverse_complement(sequence)

    strands = [("+", sequence), ("-", reverse_sequence)]
    all_orfs = []

    for strand_symbol, strand_sequence in strands:
        for frame_number in [1, 2, 3]:
            offset = frame_number - 1
            codons = get_codons(strand_sequence, frame_number)

            for start_index, stop_index in find_orfs_in_frame(codons):
                # Convert codon numbers into positions on this strand (0-based)
                start_0 = offset + start_index * 3
                end_0 = offset + (stop_index + 1) * 3

                # The ORF's DNA, including the stop codon
                orf_dna = strand_sequence[start_0:end_0]

                # Convert to 1-based positions for the user
                start_position = start_0 + 1
                end_position = end_0

                # Reverse strand: convert back to forward-strand coordinates
                if strand_symbol == "-":
                    start_position = sequence_length - end_0 + 1
                    end_position = sequence_length - start_0

                all_orfs.append({
                    "strand": strand_symbol,
                    "frame": frame_number,
                    "start": start_position,
                    "end": end_position,
                    "length": end_0 - start_0,
                    "dna": orf_dna,
                })

    return all_orfs


if __name__ == "__main__":
    test_sequences = {
        "Forward ORF": "ATGGCCAAATTTTAA",
        "Reverse ORF": "GGCTATTTCATCC",
        "Nested ATG": "ATGATGAAATAG",
        "No stop codon": "ATGGCCAAA",
    }

    for name, test_dna in test_sequences.items():
        print(f"--- {name}: {test_dna}")
        orfs = find_all_orfs(test_dna)

        if len(orfs) == 0:
            print("No ORFs found.")

        for orf in orfs:
            print(
                f"{orf['strand']}{orf['frame']} | start {orf['start']} | "
                f"end {orf['end']} | length {orf['length']} | {orf['dna']}"
            )
        print()