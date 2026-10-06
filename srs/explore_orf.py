from Bio.Seq import Seq

dna_string = "ATGGCCAAATTTTAA"

dna_seq = Seq(dna_string)

print("DNA sequence:", dna_seq)
print("Length:", len(dna_seq))

first_codon = dna_seq[0:3]
print("First codon:", first_codon)

if first_codon == "ATG":
    print("This sequence begins with a start codon!")

protein = dna_seq.translate()
print("Protein:", protein)