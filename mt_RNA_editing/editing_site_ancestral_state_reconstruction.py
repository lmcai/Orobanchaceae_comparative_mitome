from Bio import SeqIO
from ete3 import Tree
import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

#prepare phylogeny for MPR in R by naming nodes and unroot tree
t=Tree('(Aureolaria_grandiflora,((Brandisia_kwangsiensis,(((((((((Aphyllon_purpureum,(Orobanche_fasciculata,Aphyllon_fasciculatum)),(((Orobanche_ludoviciana,Orobanche_dugesii),Orobanche_cooperi),Orobanche_pinorum)),Orobanche_ramosa),((Orobanche_minor,Orobanche_crenata),Orobanche_austrohispanica)),(Mannagettae_hummelii,((Cistanche_salsa,Cistanche_deserticola),Cistanche_tubulosa))),((Epifagus_virginiana,(Conopholis_americana,Conopholis_alpina)),Kopsiopsis_strobilacea)),((((Alectra_sessiliflora,(Escobedia_crassipes,Melasma_hispidum)),Centranthera_chevalieri),((Hyobanche_sanguinea,(Aeginetia_indica,Christisonia_kwangtungensis)),Harveya_capensis)),Monochasma_japonicum)),(Lindenbergia_grandiflora,Rehmannia_glutinosa)),((Euphrasia_cuspidata,((Bartsia_inaequalis,Bellardia_latifolia),Odontites_vulgaris)),(Lathraea_squamaria,Lathraea_clandestina)))),(Pedicularis_attollens,Pedicularis_chinensis)),((Castilleja_indivisa,Castilleja_paramensis),Orthocarpus_attenuatus));')
i=1
for node in t.traverse("postorder"):
	if not node.is_leaf():
		node.name='NODE'+str(i)
		i=i+1

t.write(format=8)
'(Aureolaria_grandiflora,((Brandisia_kwangsiensis,(((((((((Aphyllon_purpureum,(Orobanche_fasciculata,Aphyllon_fasciculatum)NODE1)NODE2,(((Orobanche_ludoviciana,Orobanche_dugesii)NODE3,Orobanche_cooperi)NODE4,Orobanche_pinorum)NODE5)NODE6,Orobanche_ramosa)NODE7,((Orobanche_minor,Orobanche_crenata)NODE8,Orobanche_austrohispanica)NODE9)NODE10,(Mannagettae_hummelii,((Cistanche_salsa,Cistanche_deserticola)NODE11,Cistanche_tubulosa)NODE12)NODE13)NODE14,((Epifagus_virginiana,(Conopholis_americana,Conopholis_alpina)NODE15)NODE16,Kopsiopsis_strobilacea)NODE17)NODE18,((((Alectra_sessiliflora,(Escobedia_crassipes,Melasma_hispidum)NODE19)NODE20,Centranthera_chevalieri)NODE21,((Hyobanche_sanguinea,(Aeginetia_indica,Christisonia_kwangtungensis)NODE22)NODE23,Harveya_capensis)NODE24)NODE25,Monochasma_japonicum)NODE26)NODE27,(Lindenbergia_grandiflora,Rehmannia_glutinosa)NODE28)NODE29,((Euphrasia_cuspidata,((Bartsia_inaequalis,Bellardia_latifolia)NODE30,Odontites_vulgaris)NODE31)NODE32,(Lathraea_squamaria,Lathraea_clandestina)NODE33)NODE34)NODE35)NODE36,(Pedicularis_attollens,Pedicularis_chinensis)NODE37)NODE38,((Castilleja_indivisa,Castilleja_paramensis)NODE39,Orthocarpus_attenuatus)NODE40);'

#create binary editing status matrix
def conserved_sites_eval(i,recs):
	result=0
	sites=recs[:,i]
	if sites.count('E')>1:
		composition=list(set(sites))
		if composition==['E'] or composition==['E','-']:result=1
		if composition==['C','E'] or composition==['C','E','-']:
			if sites.count('E')>21:result=1
		if sites.count('T')>0:result=1
	return(result)

#get the list of all 43 species
recs=SeqIO.index('ccmB.deepredmt_marked.aln.fas','fasta')
all_sp=list(recs.keys())

global_converved_site=[]
genes=os.listdir('./')
genes=[i for i in genes if i.endswith('fas') and not 'exon' in i]
#initiate new seq value
new_seq_str={}
for sp in all_sp:
	new_seq_str[sp]=''

for gene in genes:
	recs=AlignIO.read(gene,'fasta')
	recs4output=SeqIO.index(gene,'fasta')
	gen_len=recs.get_alignment_length()
	#identify conserved editing sites
	converved_site=[]
	for i in range(0,gen_len):
		if conserved_sites_eval(i,recs):converved_site.append(i)
	#extract sequences of these regions
	for sp in all_sp:
		for j in converved_site:
			try:
				new_seq_str[sp]=new_seq_str[sp]+recs4output[sp].seq[j]
			except KeyError:new_seq_str[sp]=new_seq_str[sp]+'-'	
	#record these site positions
	for j in converved_site:
		global_converved_site.append(gene.split('.')[0]+'_'+str(j))

new_aln=[]
for sp in all_sp:
	seq_record = SeqRecord(Seq(new_seq_str[sp]), id=sp, description=sp)
	new_aln.append(seq_record)

new_aln = MultipleSeqAlignment(new_aln)
d=AlignIO.write(new_aln, 'conserved_editing.fasta', 'fasta')

out=open('conserved_editing.pos.txt','a')
out.write('\n'.join(global_converved_site))
out.close()


###manual examine to identify the mask regions subject to retroprcess (blocks of C->T mutations)
#This result in conserved_editing4MPR.fasta

#############################################################################################
#Direct count c->t mutations using dollo-parsimony in ete3
#Get monophyletic groups and then record the counts using node name
t=Tree('(Aureolaria_grandiflora,((Brandisia_kwangsiensis,(((((((((Aphyllon_purpureum,(Orobanche_fasciculata,Aphyllon_fasciculatum)),(((Orobanche_ludoviciana,Orobanche_dugesii),Orobanche_cooperi),Orobanche_pinorum)),Orobanche_ramosa),((Orobanche_minor,Orobanche_crenata),Orobanche_austrohispanica)),(Mannagettae_hummelii,((Cistanche_salsa,Cistanche_deserticola),Cistanche_tubulosa))),((Epifagus_virginiana,(Conopholis_americana,Conopholis_alpina)),Kopsiopsis_strobilacea)),((((Alectra_sessiliflora,(Escobedia_crassipes,Melasma_hispidum)),Centranthera_chevalieri),((Hyobanche_sanguinea,(Aeginetia_indica,Christisonia_kwangtungensis)),Harveya_capensis)),Monochasma_japonicum)),(Lindenbergia_grandiflora,Rehmannia_glutinosa)),((Euphrasia_cuspidata,((Bartsia_inaequalis,Bellardia_latifolia),Odontites_vulgaris)),(Lathraea_squamaria,Lathraea_clandestina)))),(Pedicularis_attollens,Pedicularis_chinensis)),((Castilleja_indivisa,Castilleja_paramensis),Orthocarpus_attenuatus));')
t.set_outgroup('Rehmannia_glutinosa')
for node in t.traverse("postorder"):
	if not node.is_leaf():
		node.name=0
	else:
		node.support=0
		
recs=AlignIO.read('conserved_editing4MPR.fasta','fasta')
all_sp=open('conserved_editing4MPR.fasta').readlines()
all_sp=[i[1:].strip() for i in all_sp if i.startswith('>')]
gen_len=recs.get_alignment_length()

for i in range(0,gen_len):
	#count losses of rna editing
	#editing site is not ancestral 
	if not recs[:,i][-1]=='T':
		positions = [j for j, letter in enumerate(recs[:,i]) if letter == 'T']
		if positions:
			#initiate values
			for leaf in t:leaf.add_features(editing="N")
			#map back to phylogeny
			for leaf in t:
				if leaf.name in [all_sp[j] for j in positions]:leaf.add_features(editing="T")
			for node in t.get_monophyletic(values="T", target_attr="editing"):
				#internal branch
				if not node.is_leaf():node.name=node.name+1
				#terminal branch
				else:node.support=node.support+1

#terminal branch
for leaf in t:print(leaf.name,leaf.support)

Rehmannia_glutinosa 0.0
Lindenbergia_grandiflora 13.0
Aphyllon_purpureum 2.0
Orobanche_fasciculata 0.0
Aphyllon_fasciculatum 0.0
Orobanche_ludoviciana 0.0
Orobanche_dugesii 0.0
Orobanche_cooperi 1.0
Orobanche_pinorum 0.0
Orobanche_ramosa 8.0
Orobanche_minor 0.0
Orobanche_crenata 0.0
Orobanche_austrohispanica 3.0
Mannagettae_hummelii 7.0
Cistanche_salsa 1.0
Cistanche_deserticola 0.0
Cistanche_tubulosa 9.0
Epifagus_virginiana 2.0
Conopholis_americana 6.0
Conopholis_alpina 2.0
Kopsiopsis_strobilacea 3.0
Alectra_sessiliflora 2.0
Escobedia_crassipes 0.0
Melasma_hispidum 0.0
Centranthera_chevalieri 2.0
Hyobanche_sanguinea 2.0
Aeginetia_indica 3.0
Christisonia_kwangtungensis 2.0
Harveya_capensis 4.0
Monochasma_japonicum 7.0
Euphrasia_cuspidata 1.0
Bartsia_inaequalis 0.0
Bellardia_latifolia 5.0
Odontites_vulgaris 0.0
Lathraea_squamaria 1.0
Lathraea_clandestina 1.0
Brandisia_kwangsiensis 1.0
Pedicularis_attollens 0.0
Pedicularis_chinensis 1.0
Aureolaria_grandiflora 1.0
Castilleja_indivisa 0.0
Castilleja_paramensis 2.0
Orthocarpus_attenuatus 1.0
				
#internal counts
t.write(format=8)
'(Rehmannia_glutinosa,(Lindenbergia_grandiflora,((((((((Aphyllon_purpureum,(Orobanche_fasciculata,Aphyllon_fasciculatum)2)1,(((Orobanche_ludoviciana,Orobanche_dugesii)0,Orobanche_cooperi)0,Orobanche_pinorum)2)3,Orobanche_ramosa)5,((Orobanche_minor,Orobanche_crenata)1,Orobanche_austrohispanica)15)3,(Mannagettae_hummelii,((Cistanche_salsa,Cistanche_deserticola)0,Cistanche_tubulosa)4)4)0,((Epifagus_virginiana,(Conopholis_americana,Conopholis_alpina)2)4,Kopsiopsis_strobilacea)1)0,((((Alectra_sessiliflora,(Escobedia_crassipes,Melasma_hispidum)0)0,Centranthera_chevalieri)0,((Hyobanche_sanguinea,(Aeginetia_indica,Christisonia_kwangtungensis)10)0,Harveya_capensis)2)6,Monochasma_japonicum)2)0,(((Euphrasia_cuspidata,((Bartsia_inaequalis,Bellardia_latifolia)0,Odontites_vulgaris)1)10,(Lathraea_squamaria,Lathraea_clandestina)0)6,(Brandisia_kwangsiensis,((Pedicularis_attollens,Pedicularis_chinensis)3,(Aureolaria_grandiflora,((Castilleja_indivisa,Castilleja_paramensis)0,Orthocarpus_attenuatus)1)3)6)0)1)0)0);'
	

#convert c->t mutations to 1 and everything else is 0 (including indels, other mutations)
recs=AlignIO.read('conserved_editing4MPR.fasta','fasta')
all_sp=open('conserved_editing4MPR.fasta').readlines()
all_sp=[i[1:].strip() for i in all_sp if i.startswith('>')]
gen_len=recs.get_alignment_length()

Tmut={}
for sp in all_sp:
	Tmut[sp]=[]


for i in range(0,gen_len):	
	#editing site is not ancestral 
	if not recs[:,i][-1]=='T':
		positions = [j for j, letter in enumerate(recs[:,i]) if letter == 'T']
		if positions:
			for k in Tmut.keys():
				if k in [all_sp[j] for j in positions]:Tmut[k].append(1)
				else:Tmut[k].append(0)

out=open('MPR.01matrix.tsv','a')
for sp in all_sp:
	d=out.write(sp+'\t'+'\t'.join([str(i) for i in Tmut[sp]])+'\n')

out.close()

###execute MPR.R to obtain MPR_lower.csv and MPR_upper.csv
#To count changes based on ancestral state reconstruction
from ete3 import Tree
t=Tree('MPR_output.tre',format=8)
all_states=open('MPR_lower.csv').readlines()
number_sites=len(all_states[0].split(','))
Tmut={}
for i in range(1,number_sites):
	states={}
	for l in all_states:states[l.split(',')[0]]=l.split(',')[i]
	for node in t.traverse():
		if not node.name=='Root' and not node.is_leaf():
			if states[node.get_ancestors()[0].name]!=states[node.name] and states[node.name]=='1':
				try:Tmut[node.name]=Tmut[node.name]+1
				except KeyError:Tmut[node.name]=1

for node in t.traverse():
	if node.name in list(Tmut.keys()):
		node.name=Tmut[node.name]

t.write(format=8)

'(((((((Aureolaria_grandiflora,((Castilleja_indivisa,Castilleja_paramensis)N39,Orthocarpus_attenuatus)4)6,(Pedicularis_attollens,Pedicularis_chinensis)8)N36,Brandisia_kwangsiensis)1,((Euphrasia_cuspidata,((Bartsia_inaequalis,Bellardia_latifolia)N30,Odontites_vulgaris)1)2,(Lathraea_squamaria,Lathraea_clandestina)N33)6)1,(((((((Aphyllon_purpureum,(Orobanche_fasciculata,Aphyllon_fasciculatum)2)1,(((Orobanche_ludoviciana,Orobanche_dugesii)N3,Orobanche_cooperi)N4,Orobanche_pinorum)2)3,Orobanche_ramosa)5,((Orobanche_minor,Orobanche_crenata)1,Orobanche_austrohispanica)15)2,(Mannagettae_hummelii,((Cistanche_salsa,Cistanche_deserticola)N11,Cistanche_tubulosa)4)4)N14,((Epifagus_virginiana,(Conopholis_americana,Conopholis_alpina)2)3,Kopsiopsis_strobilacea)N17)N18,((((Alectra_sessiliflora,(Escobedia_crassipes,Melasma_hispidum)N19)N20,Centranthera_chevalieri)N21,((Hyobanche_sanguinea,(Aeginetia_indica,Christisonia_kwangtungensis)10)N23,Harveya_capensis)2)6,Monochasma_japonicum)1)1)N28,Lindenbergia_grandiflora)1,Rehmannia_glutinosa);'