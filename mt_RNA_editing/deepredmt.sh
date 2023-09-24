deepredmt test.fas > seqs.pred
cut -f3,5 seqs.pred | sed 's/C/0/g' | sed 's/E/1/g' > seqs.parsed_pred

#Then execute across genes
deepredmt atp1.aln.fas >atp1.deepredmt.pred

#Then summarize
awk '{if ($5 >0.9) print FILENAME,$1,$5}' *.deepredmt.pred >deepredmt.all.pred  

#Get sum stat in python
!/usr/bin/python
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
		except KeyError:editing_site_num.append(0)
	editing_site_num=[str(i) for i in editing_site_num]
	out.write(k+'\t'+'\t'.join(editing_site_num)+'\n')
