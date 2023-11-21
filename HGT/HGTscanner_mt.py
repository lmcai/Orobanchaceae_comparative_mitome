from Bio import SeqIO
import os,argparse
import datetime
from ete3 import Tree
import ete3
import pybedtools
from numpy import median


parser = argparse.ArgumentParser(description='HGTscanner_mt is a utility wrappepr to identify HGT blocks in organellar (mostly mitochondrial) genomes.')
parser.add_argument('-q', metavar='query', help='fasta file of the target mitome', required=True)
parser.add_argument('-ref', metavar='reference', help='one fasta file containing all custom references of both close relatives and potential HGT donor. This will be combined with the NCBI Viridiplantae mito database.')
parser.add_argument('-o', metavar='output', help='output prefix', required=True)
parser.add_argument('-f', metavar='family', help='family of the query for HGT classification', required=True)
parser.add_argument('-b', metavar='mask', help='bed file for regions to be masked')

#pass argument values to variables
args = parser.parse_args()
query=args.q
if args.ref:
	reference=args.ref
sp=args.o
fam=args.family
if args.b:
	mask_bed=args.b


###############
#MASK gene/MTPT regions
def mask_fasta_with_bed(fasta_file, bed_file, output_file):
	sequences = SeqIO.index(fasta_file, "fasta")
	bed = open(bed_file).readlines()
	recs={}
	for k in sequences.keys():
		recs[k]=str(sequences[k].seq)
	for l in bed:
   		sequence_id = l.split()[0]
   		start = int(l.split()[1])
   		end = int(l.split()[2])
   		sequence = recs[sequence_id]
   		masked_sequence = sequence[:(start-1)] + "N" * (end - start + 1) + sequence[end:]
   		recs[sequence_id] = masked_sequence
	out=open(output_file,'w')
	for k in recs.keys():
		d=out.write('>'+k+'\n'+recs[k]+'\n')



if args.b:
	print(str(datetime.datetime.now())+'\tMasking query sequence '+query+' using the bed file: '+args.b)
	mask_fasta_with_bed(query, mask_bed, sp+'.mask.fas')


###########
#BLAST
if args.ref:
	#add custom references to database
	print(str(datetime.datetime.now())+'\tAdded custom reference '+reference+' to the NCBI Viridiplantae mitochondrial database')
	S='cat '+reference+' Viridiplantae_mt.fasta >mt_db.fas'
	os.system(S)
else:
	S='cp Viridiplantae_mt.fasta mt_db.fas'
	os.system(S)


S='makeblastdb -in mt_db.fas -out mt -dbtype nucl >/dev/null'
os.system(S)


if args.b:S='blastn -task dc-megablast -query '+sp+'.mask.fas -db mt -outfmt 6 -evalue 1e-20 >'+sp+'.raw.blast'
else:S='blastn -task dc-megablast -query '+query+' -db mt -outfmt 6 -evalue 1e-20 >'+sp+'.raw.blast'
os.system(S)
print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)


####################################
#Add taxonomic information to each hit at species and family level
x=open(sp+'.raw.blast').readlines()
out=open(sp+'.taxon.blast','w')
species={}
family={}
y=open('Viridiplantae_mt.taxonomy').readlines()
for l in y:
	species[l.split('\t')[0]]=l.split('\t')[1]
	family[l.split('\t')[0]]=l.split('\t')[2].strip()


for l in x:
	try:
		d=out.write(l.strip()+'\t'+species[l.split()[1]]+'\t'+family[l.split()[1]]+'\n')
	except KeyError:
		if l.split()[1].startswith('Oro'):d=out.write(l.strip()+'\t'+l.split()[1]+'\tOrobanchaceae\n')
		else:d=out.write(l.strip()+'\t'+l.split()[1]+'\tNA\n')


out.close()

###########
#sort blast results, give each row a unique id, only hits <20k are printed
S="cat "+sp+".taxon.blast | sort -k1,1 -k7,7n -k11,11n | awk -v OFS='\\t' '{if ($8-$7 <= 20000) print $1, $7, $8, $2, $9, $10, $13, $14, NR}' > "+sp+".taxon_sorted.bed"
os.system(S)

####################################
#define potential HGT blocks

#classify these hits based on source, here, alignments from Orobanchaceae, Lamiaceae, Oleaceae, and Scrophulariaceae will be regarded as close relatives
#This is necessary because they create long synteny blocks that may be broken by HGT in the middle

x=open(sp+".taxon_sorted.bed").readlines()
otherfam=[]
samefam=[]
for l in x:
	if not l.split()[7] in ['Orobanchaceae','Lamiaceae','Scrophulariaceae','Oleaceae']:
		otherfam.append(l)
	else:
		samefam.append(l)

otherfam_merged=pybedtools.BedTool(''.join(otherfam), from_string=True).merge(c=9,o='collapse')
samefam_bed=pybedtools.BedTool(''.join(samefam), from_string=True)
out=open(sp+'.merged.bed','w')
d=out.write(str(otherfam_merged))

print(str(datetime.datetime.now())+'\tFound '+str(len(otherfam_merged))+' homologous genetic blocks for further exaination')

###############everything below in this section is abandoned
def id2bed(ids,bed_file):
	id_dict={}
	for l in bed_file:
		id_dict[l.split()[-1]]=l
	filtered_bed=[id_dict[j] for j in ids]
	return(filtered_bed)

def bedOK2go(bed_txt,all_aln):
	ids=bed_txt.split('\t')[-1].strip()
	aln=id2bed(ids.split(','),all_aln)
	median_len=median([int(j.split()[2])-int(j.split()[1]) for j in aln])
	bed_start=min([int(j.split()[1]) for j in aln])
	bed_end=max([int(j.split()[2]) for j in aln])
	bed_len=bed_end - bed_start
	#get median blast hit length, if <500 bp or largely overlap with the bed size, output
	if median_len>0.5*bed_len or bed_len<500:return(1)
	else:return(0)

def break_up_range(start, end, step):
    ranges = []
    for i in range(start, end + 1, step):
        ranges.append((i, min(i + step - 1, end)))
    return ranges


#This function is abandoned since I did not find a good way to cluster individual hits to groups with minor differences on their boundries
def break_large_beds(bed_file):
	#I. remove blast hits from the same family to build new merged bed file (hopefully they are smaller). These hits will be added back don't worry
	filtered=[]
	addback=[]
	for j in bed_file:
		if not (j.split('\t')[7]==fam or j.split('\t')[3].startswith('Oro')):
			filtered.append(j)
		else:
			addback.append(j)
	filtered_bed = pybedtools.BedTool(''.join(filtered), from_string=True).merge(c=9,o='collapse')
	#if merged_bed meet the criteria, you can define interactions with '+' and '-'!!
	filtered_bed_str=str(filtered_bed).strip()
	OK=1
	for l in filtered_bed_str.split('\n'):
		if not bedOK2go(l,y):OK=0
	if OK:
		#removing same family aln is sufficient to create acceptable homologs
		addback_ids=[j.split()[-1] for j in addback]
		#add back aln ids to each smaller bed region 
		filtered_bed_str=filtered_bed_str.strip()
		output_txt=''
		for l in filtered_bed_str.split('\n'):
			output_txt=output_txt+l+','+','.join(addback_ids)+'\n'
		return(output_txt)
	else:
		#II. break up blocks based on coverage patterns
		return('')
		#full_bed=str(pybedtools.BedTool(''.join(bed_file), from_string=True).merge())
		#broken_range=break_up_range(int(full_bed.split()[1]),int(full_bed.split()[2]),5)
		#broken_bed_txt=[full_bed.split()[0]+'\t'+str(i[0])+'\t'+str(i[1]) for i in broken_range]
		#broken_bed=pybedtools.BedTool('\n'.join(broken_bed_txt), from_string=True)
		#all_bed=pybedtools.BedTool(''.join(bed_file), from_string=True)
		#bed_cov=broken_bed.coverage(all_bed)
		
		#filtered_bed_len={}
		#filtered_bed_len_list=[]
		#for j in filtered:
		#	try:
		#		filtered_bed_len[str(int(j.split()[2])-int(j.split()[1]))].append(j)
		#	except KeyError:filtered_bed_len[str(int(j.split()[2])-int(j.split()[1]))]=[j]
		#	filtered_bed_len_list.append(int(j.split()[2])-int(j.split()[1]))
		#while True:
		#	cur_max=max(filtered_bed_len_list)
		#	cur_max_bed=filtered_bed_len[cur_max]
		#	for k in cur_max_bed:
		#		filtered.remove(k)
		#		addback.append(k)
		#		filtered_bed_len_list.remove(cur_max)
		#		filtered_bed = pybedtools.BedTool(''.join(filtered), from_string=True).merge()
		#		if bedOK2go(str(filtered_bed)):
					#save
		#			break

####################################
#extract sequences for each block

q_recs=SeqIO.index(query,'fasta')
ref_recs=SeqIO.index('mt_db.fas', 'fasta')

i=1
for hit in otherfam_merged:
	#gather overlapping alignment from both other families and close relatives
	ids=hit.fields[3]
	raw_beds=id2bed(ids.split(','),y)
	out=open(sp+'.'+str(i)+'.fas','w')
	for l in raw_beds:
		d=out.write('>'+'|'.join(l.split()[3:8])+'\n')
		d=out.write(str(ref_recs[l.split()[3]].seq[(int(l.split()[4])-1):int(l.split()[5])])+'\n')
		
	#write query
	d=SeqIO.write(q_recs[hit.chrom][(hit.start-1):hit.end],out,'fasta')
	#add hits overlap with close relatives
	
	i=i+1
	out.close()

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