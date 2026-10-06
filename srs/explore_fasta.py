from pathlib import Path
from Bio import SeqIO

project_folder = Path(__file__).resolve().parent.parent
fasta_path = project_folder / "data" / "sample.fasta"


for record in SeqIO.parse(fasta_path, "fasta"):
    print("ID:", record.id)
    print("Description:", record.description)
    print("Length:", len(record.seq))
    print("Sequence:", record.seq)
    print()