from ete3 import Tree
import ete3
import os,sys

close_relative={'Orobanchaceae','Lamiaceae','Bignoniaceae','Phrymaceae','Acanthaceae','Oleaceae','Scrophulariaceae','Verbenaceae','Plantaginaceae','Gesneriaceae','Lentibulariaceae'}
id2sp={'Ain':"Aeginetia_indica",
"Ase":"Alectra_sessiflora",
"Afa":"Aphyllon_faciculatum",
"Apr":"Aphyllon_purpureum",
"Agr":"Aureolaria_grandiflora",
"Bin":"Bartsia_inaequalis",
"Bla":"Bellardia_latifolia",
"Bkw":"Brandisia_kwangsiensis",
"Cin":"Castilleja_indivisa",
"Cpa":"Castilleja_paramensis",
"Cch":"Centranthera_chevalieri",
"Ckw":"Christisonia_kwangtungensis",
"Cde":"Cistanche_deserticola",
"Csa":"Cistanche_salsa",
"Ctu":"Cistanche_tubulosa",
"Cal":"Conopholis_alpina",
"Cam":"Conopholis_americana",
"Evi":"Epifagus_virginiana",
"Ecr":"Escobedia_crassipes",
"Ecu":"Euphrasia_cuspidata",
"Hca":"Harveya_capensis",
"Hsa":"Hyobanche_sanguinea",
"Kst":"Kopsiopsis_strobilacea",
"Lgr":"Lindenbergia_grandiflora",
"Mhu":"Mannagettaea_hummelii",
"Mhi":"Melasma_hispidum",
"Mja":"Monochasma_japonicum",
"Ose":"Odontites_serotina",
"Oau":"Orobanche_austrohispanica",
"Oco":"Orobanche_cooperi",
"Ocr":"Orobanche_crenata",
"Odu":"Orobanche_dugesii",
"Ofa":"Orobanche_faciculata",
"Ofo":"Orobanche_foetida",
"Olu":"Orobanche_ludoviciana",
"Omi":"Orobanche_minor",
"Opi":"Orobanche_pinorum",
"Ora":"Orobanche_ramosa",
"Oat":"Orthocarpus_attenuatus",
"Pat":"Pedicularis_attollens",
"Pch":"Pedicularis_chinensis",
"Rgl":"Rehmannia_glutinosa",
}

#take a list of fasta headers or tree tips, return species names
def ID2sp(lst,target_sp):
	outputsp=[]
	for j in lst:
		sp=j.split('|')[-1]
		if sp==target_sp:
			sp=id2sp[sp.split('_')[0]]
		elif sp.startswith('Oro') and not sp.startswith('Orobanche'):
			sp=sp[3:]
			sp=id2sp[sp]
		outputsp.append(sp)
	return(list(set(outputsp)))


# Function to find the largest continuous block in a list around a specific element based on its name
def find_max_block_around_element(lst, taxonomy_name, element_index, max_diff_count):
	diff_count = 0
	current_start = element_index
	current_end = element_index
	first_ingroup=element_index
	last_ingroup=element_index
	ingroup_sp=[]
	outgroup_sp=[]
	alien_fam='NA'
	#extend on the left side
	for i in range(current_start,0,-1):
		tip, taxonomy = lst[i]
		if taxonomy == taxonomy_name:  # When same as the target taxon rank
			current_start = max(0, current_start - 1)  # Start searching from one element before
			diff_count=0
			first_ingroup=i
			ingroup_sp.append(tip)
		else:
			if not taxonomy == alien_fam:#a new, non-orobanchaceae family
				if diff_count< max_diff_count:current_start = max(0, current_start - 1)
				else:break
				diff_count += 1
				alien_fam=taxonomy
				if not taxonomy in close_relative:outgroup_sp.append(tip)
			else:#a same, non-orobanchaceae family
				current_start = max(0, current_start - 1)
				if not taxonomy in close_relative:outgroup_sp.append(tip)
	#extend on the right side
	for i in range(current_end,len(lst)):
		tip, taxonomy = lst[i]
		if taxonomy == taxonomy_name:  # When same as the target taxon rank
			current_end = min(len(lst), current_end+1)  # Start searching from one element before
			diff_count=0
			last_ingroup=i
			ingroup_sp.append(tip)
		else:
			if not taxonomy == alien_fam:#a new, non-orobanchaceae family
				if diff_count< max_diff_count:current_end= min(len(lst), current_end+1)
				else:break
				diff_count += 1
				alien_fam=taxonomy
				if not taxonomy in close_relative:outgroup_sp.append(tip)
			else:#a same, non-orobanchaceae family
				current_start = max(0, current_start - 1)
				if not taxonomy in close_relative:outgroup_sp.append(tip)
	return(first_ingroup, last_ingroup, ingroup_sp, outgroup_sp)

def VGTFromTipOrder(tree,target_sp):
	tips=[node.name for node in tree]
	allfamilies=[]
	orders=[]
	for j in tips:
		if j==target_sp:
			allfamilies.append((j,'Orobanchaceae'))
			orders.append('Lamiales')
		else:
			allfamilies.append((j,j.split('|')[0]))
			if j.split('|')[0] in close_relative:orders.append('Lamiales')
			else:orders.append('NOT')
	max_different_names = 2  # Maximum number of different names allowed
	element_index = next(index for index, (name, _) in enumerate(allfamilies) if name == target_sp)
	oro_start,oro_end,oro_tip=find_max_block_around_element(allfamilies, 'Orobanchaceae', element_index, max_different_names)
	oro_genera=[j.split('_')[0] for j in ID2sp(oro_tip,target_sp)]
	oro_genera=set(oro_genera)
	#if this orobanchaceae cluster contains more than five genera, this is most likely a VGT
	if len(oro_genera)>4:
		return(1)
	#check if this cluster is nested within a larger Lamiales cluster
	else:
		output=0
		if oro_start!=0:
			if orders[oro_start-1]=='Lamiales':output=1
		if oro_end!=(len(tips)-1):
			if orders[oro_end+1]=='Lamiales':output=1
		return(output)
	
def main(target,trfile):
	t=Tree(target+'_HGTscanner_supporting_files/'+trfile)
	allsp=[node.name for node in t]
	target_tip=[j for j in allsp if j.startswith(target)]
	target_sp=target_tip[0]
	allfamilies=[]
	for j in allsp:
		if j==target_sp:
			allfamilies.append((j,'Orobanchaceae'))
		else:
			allfamilies.append((j,j.split('|')[0]))
	max_different_names = 2  # Maximum number of different names allowed
	element_index = next(index for index, (name, _) in enumerate(allfamilies) if name == target_sp)
	oro_start,oro_end,oro_tip,out_tip=find_max_block_around_element(allfamilies, 'Orobanchaceae', element_index, max_different_names)
	d=out.write(trfile+'\t'+';'.join(out_tip)+'\n')

trees=os.listdir(sys.argv[1]+'_HGTscanner_supporting_files/')
trees=[j for j in trees if j.endswith('treefile')]
out=open('alien_nested_in_Orobanchaceae.tsv','a')
for tr in trees:
	try:main(sys.argv[1],tr)
	except ete3.parser.newick.NewickError:print(tr)
out.close()