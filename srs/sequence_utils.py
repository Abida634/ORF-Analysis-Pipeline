from Bio import SeqIO
VALID_BASES = set("ATGC")


def clean_sequence(raw_sequence):
    no_whitespace = "".join(raw_sequence.split())
    return no_whitespace.upper()


def validate_sequence(sequence):
    if len(sequence) == 0:
        raise ValueError("Sequence is empty.")

    invalid_characters = set(sequence) - VALID_BASES

    if len(invalid_characters) > 0:
        invalid_list = ", ".join(sorted(invalid_characters))
        raise ValueError(f"Invalid characters found: {invalid_list}")


def prepare_dna(raw_sequence):
    cleaned = clean_sequence(raw_sequence)
    validate_sequence(cleaned)
    return cleaned

def read_fasta(file_path):

    sequences = []

    for record in SeqIO.parse(file_path, "fasta"):
        raw_sequence = str(record.seq)
        try:
            clean_dna = prepare_dna(raw_sequence)
        except ValueError as error:

            raise ValueError(f"Record '{record.id}': {error}")
        sequences.append((record.id, clean_dna))

    if len(sequences) == 0:
        raise ValueError("No FASTA records found. Is this a valid FASTA file?")

    return sequences



if __name__ == "__main__":
    test_inputs = [
        "  atg gcc\naaa ttt\ttaa \n",
        "ATGXC1T",
        "   ",
    ]

    for test in test_inputs:
        try:
            result = prepare_dna(test)
            print("OK:", result)
        except ValueError as error:
            print("Rejected:", error)