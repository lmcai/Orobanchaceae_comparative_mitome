from Bio import SeqIO

ath=open('ATH_ref.fas').readlines()
ath=[i[1:].split('.')[0] for i in ath if i.startswith('>')]

import os
ath=os.listdir('.')
ath=[i.split('.')[0] for i in ath if i.endswith('.aln.fas')]

lib1={}
lib2={}
b1=open('Phtheirospermum.1.blast').readlines()
cur='NA'
for l in b1:
	if l.split()[0]!=cur:
		lib1[l.split()[0]]=l.split()[1]
		cur=l.split()[0]



b2=open('Phtheirospermum.2.blast').readlines()
cur='NA'
for l in b2:
	if l.split()[0]!=cur:
		ath_id=l.split()[1]
		ath_id=ath_id.split('.')[0]
		lib2[l.split()[0]]=ath_id
		cur=l.split()[0]




recs=SeqIO.index('CDS/Phtheirospermum_trinity.Trinity.fasta','fasta')
for i in ath:
	try:
		if lib2[lib1[i]]==i:
			out=open(i+'.fas','a')
			d=SeqIO.write(recs[lib1[i]],out,'fasta')
			out.close()
		else:print('Not best:' +i)
	except KeyError:print('NO hit: '+ i)



recs=SeqIO.parse('ATH_ref.fas','fasta')
for rec in recs:
	out=open(rec.id.split('_')[0]+'.fas','a')
	d=SeqIO.write(rec,out,'fasta')
	out.close()