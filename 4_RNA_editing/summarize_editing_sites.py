import os, re
from Bio import SeqIO

#summarize number of RNA editing sites per species
x=open('deepredmt.all.pred').readlines()
a={}
for l in x:
	sp=l.split()[1]
	sp=sp.split('!')[0]
	gene=l.split('.')[0]
	try:
		a[sp][gene]=a[sp][gene]+1
	except KeyError:
		try:a[sp][gene]=1
		except KeyError:
			a[sp]={}
			a[sp][gene]=1

genes=[l.split('.')[0] for l in x]
genes=list(set(genes))
out=open('deepredmt.sum_stat.tsv','a')
out.write('\t'.join(['sp']+genes)+'\n')
for k in a.keys():
	editing_site_num=[]
	for j in genes:
		try:editing_site_num.append(a[k][j])
		except KeyError:editing_site_num.append('NA')
	editing_site_num=[str(i) for i in editing_site_num]
	out.write(k+'\t'+'\t'.join(editing_site_num)+'\n')

#summarize the codon position of RNA editing sites
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
