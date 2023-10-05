import os, re
from Bio import SeqIO
files=os.listdir('./')
files=[i for i in files if i.endswith('marked.aln.fas') and not i.startswith('r')]

ed1=0
ed2=0
ed3=0
for file in files:
	recs=SeqIO.index(file,'fasta')
	if len(recs['Orobanche_ludoviciana'].seq)%3==0:
		#codon is correct
		for k in recs.keys():
			positions = [match.start() for match in re.finditer('E', str(recs[k].seq))]
			for j in positions:
				if j%3==1:ed2=ed2+1
				elif j%3==2:ed3=ed3+1
				elif j%3==0:ed1=ed1+1
	else:
		print(file)


>>> ed1
4742
>>> ed2
10389
>>> ed3
1076
