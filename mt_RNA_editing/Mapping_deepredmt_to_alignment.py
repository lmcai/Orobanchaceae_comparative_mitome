import os
from Bio import SeqIO
files=os.listdir('.')
files=[i for i in files if i.endswith('pred')]

for file in files:
	x=open(file).readlines()
	recs=SeqIO.index(file.split('.pred')[0]+'.fasta','fasta')
	ed_sites={}
	for l in x:
		if float(l.split()[4])>0.9:
			sp=l.split()[0]
			try:ed_sites[sp.split('!')[0]].append(sp.split('!')[1])
			except KeyError:ed_sites[sp.split('!')[0]]=[sp.split('!')[1]]
	out=open(file.split('.')[0]+'.deepredmt_marked.fas','a')
	for k in recs.keys():
		if k in ed_sites.keys():
			newseq=str(recs[k].seq)
			for i in ed_sites[k]:
				if recs[k].seq[int(i)-1]=='C':
					newseq=newseq[:int(i)-1]+'E'+newseq[int(i):]
			d=out.write('>'+k+'\n'+newseq+'\n')
		else:
			d=SeqIO.write(recs[k],out,'fasta')
