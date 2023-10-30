from Bio import SeqIO
import os,argparse
import datetime
from ete3 import Tree
import ete3

parser = argparse.ArgumentParser(description='HGTfinder is a utility wrappepr to identify HGT blocks in organellar (mostly mitochondrial) genomes.')
parser.add_argument('-q', metavar='query', help='fasta file of the target mitome', required=True)
parser.add_argument('-ref', metavar='reference', help='one fasta file containing all references including close relatives and potential HGT donor', required=True)
parser.add_argument('-o', metavar='output', help='output prefix')

args = parser.parse_args()


query=args.q
reference=args.ref
sp=args.o

###########
#BLAST

S='makeblastdb -in '+query+' -out '+sp+' -dbtype nucl >/dev/null'
os.system(S)
S='blastn -task dc-megablast -query '+reference+' -db '+sp+' -outfmt 6 -evalue 1e-40 >'+sp+'.mtpt.blast'
os.system(S)
print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)
###########
#sort blast results and five each row an uniq id
S="awk -v OFS='\\t' '{if ($9 <= $10) print $2, $9, $10, $1, $7, $8, $3, S11, $12; else print $2, $10, $9, $1, $7, $8, $3, S11, $12}' "+sp+".mtpt.blast| sort -k1,1 -k2,2n -k4,4n | awk 'BEGIN{FS=OFS='\\t'} {print $0, NR}' > "+sp+".mtpt.bed"
#print(S)
os.system(S)

###########
#define potential HGT blocks
S="bedtools merge -i "+sp+".mtpt.bed -c 9 -o collapse >"+sp+'.temp.bed'
os.system(S)

###########
#extract sequences for each block
loci=open(sp+'.temp.bed').readlines()
hits=open(sp+'.mtpt.bed').readlines()
q_recs=SeqIO.index(query,'fasta')
ref_recs=SeqIO.index(reference, 'fasta')

seq_loc={}
for l in hits:
	seq_loc[l.split()[-1]]=l

order=1
for l in loci:
	out=open(sp+'.temp.'+str(order)+'.fas','w')
	d=SeqIO.write(q_recs[l.split()[0]][(int(l.split()[1])-1):int(l.split()[2])],out,'fasta')
	out.close()
	out=open(sp+'.temp.'+str(order)+'.fas','a')
	bed_ids=l.split()[3].strip()
	bed_ids=bed_ids.split(',')
	for i in bed_ids:
		line=seq_loc[i]
		id=line.split()[3]
		start=int(line.split()[4])-1
		end=int(line.split()[5])
		d=out.write('>'+id+'_'+str(start)+'_'+str(end)+'\n')
		d=out.write(str(ref_recs[id].seq[start:end])+'\n')
	out.close()
	order=order+1

print(str(datetime.datetime.now())+'\tExatracted potential HGT sequences for '+sp)

#############
#alignment and phylogenetic reconstruction
print(str(datetime.datetime.now())+'\tStart alignment and phylogenetic reconstruction with mafft and iqtree for '+str(order-1)+' regions. May take a while...')

for i in range(1,order):
	S="mafft --quiet --adjustdirection "+sp+".temp."+str(i)+".fas | sed 's/_R_//g' > "+sp+".gt."+str(i)+".aln.fas"
	os.system(S)
	b=SeqIO.index(sp+".gt."+str(i)+".aln.fas",'fasta')
	q=loci[i-1].split()[0]
	if len(b[q].seq)<10000:
		S="nohup iqtree -B 1000 -T 4 --quiet -m GTR+F -redo -s "+sp+".gt."+str(i)+".aln.fas >/dev/null 2>&1"
		os.system(S)
		print(str(datetime.datetime.now())+'\tLoci #'+str(i))
	else:print(str(datetime.datetime.now())+'\tLoci #'+str(i)+' is longer than 10kb. Skip tree building. Check manually.')

os.system('rm '+sp+'*.bionj')
os.system('rm '+sp+'*.gz')
os.system('rm '+sp+'*.log')
os.system('rm '+sp+'*.iqtree')
os.system('rm '+sp+'*.mldist')
os.system('rm '+sp+'*.phy')
os.system('rm '+sp+'*.contree')
os.system('rm '+sp+'*.nex')

##############
#Evaluate the source of the region and output summary file
out=open(sp+'.hgt.sum.tsv','w')
out.write('ID\tTarget_scaffold\tStart\tEnd\tPhylo_source\tBlast_hit_ID\n')
out.close()
out=open(sp+'.hgt.sum.tsv','a')

for i in range(1,order):
	q=loci[i-1].split()[0]
	outgroup=[]
	try:
		t=Tree(sp+'.gt.'+str(i)+'.aln.fas.treefile')
		ancestor=t.get_midpoint_outgroup()
		t.set_outgroup(ancestor)
		q_branch=t&q
		if not q_branch.get_ancestors()[0].is_root():
			sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
		else:
			t=Tree(sp+'.gt.'+str(i)+'.aln.fas.treefile')
			outgroup=[leaf.name for leaf in t if leaf.name.startswith('Sorghum')]
			if len(outgroup)>0:
				t.set_outgroup(outgroup[0])
				q_branch=t&q
				sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
			else:
				outgroup=[leaf.name for leaf in t if leaf.name.startswith('Rehm')]
				if outgroup:
					t.set_outgroup(t&outgroup[0])
					q_branch=t&q
					sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
				else:
					#no sorghum no rehmannia
					sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
		out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
	except ete3.parser.newick.NewickError:
		sisters=open(sp+'.temp.'+str(i)+".fas").readlines()
		sisters=[i[1:].strip() for i in sisters if (i.startswith('>')) and (not i[1:].strip()==q)]
		out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+'ALL HITS: '+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
		
out.close()
os.system('rm '+sp+'.temp.*.fas')
os.system('rm '+sp+'.temp.bed')
os.system('rm '+sp+'.n*')

if not os.path.isdir('HGTscanner_supporting_files'):os.mkdir('HGTscanner_supporting_files')
os.system('mv '+sp+'.*.aln.fas* HGTscanner_supporting_files')
print(str(datetime.datetime.now())+'\tCompleted evaluation of HGT source. See summary file in '+sp+'.hgt.sum.tsv')