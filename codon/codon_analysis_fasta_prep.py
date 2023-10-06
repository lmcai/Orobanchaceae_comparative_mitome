from Bio import SeqIO
import os
files=os.listdir('./')
files=[i for i in files if i.endswith('.fasta') and not i.startswith('r')]

for file in files:
	recs=SeqIO.parse(file,'fasta')
	for rec in recs:
		out=open(rec.id+'.codon.fas','a')
		d=out.write('>'+file.split('.')[0]+'')
		d=out.write('\n'+str(rec.seq)+'\n')
		out.close()

exit()


#filter seqeunce for frameshifts
from Bio import SeqIO

import os
files=os.listdir('./')
files=[i for i in files if i.endswith('.codon.fas')]

def count_continuous_dashes(sequence):
    counts = []
    current_count = 0
    for char in sequence:
        if char == '-':
            current_count += 1
        elif current_count > 0:
            counts.append(current_count)
            current_count = 0
    if current_count > 0:
        counts.append(current_count)
    return counts

for file in files:
	print('\n\n'+file)
	for record in SeqIO.parse(file, "fasta"):
		sequence = str(record.seq)
		counts = count_continuous_dashes(sequence)
		print(f"{record.id}: {', '.join(map(str, counts))}")
