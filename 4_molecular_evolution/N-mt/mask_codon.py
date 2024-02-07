from Bio import SeqIO
from Bio.Seq import Seq
import sys

def mask_invalid_codons(sequence):
    masked_sequence = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        if len(codon) == 3:
            try:
                codon_seq = Seq(codon)
                aa = codon_seq.translate()
                masked_sequence += codon if (aa != "*" and aa != "X") else "---"
            except Exception as e:
                masked_sequence += "NNN"
        else:
            masked_sequence += codon
    return masked_sequence

# Read the input FASTA file
input_file = sys.argv[1]
output_file = open(sys.argv[2],'w')

for record in SeqIO.parse(input_file, "fasta"):
	if not record.id.startswith('Reh'):
		masked_sequence = mask_invalid_codons(str(record.seq))
		record.seq = masked_sequence
		output_file.write('>'+record.id+'\n'+masked_sequence+'\n')
	

output_file.close()

