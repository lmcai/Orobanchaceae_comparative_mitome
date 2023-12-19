from ete3 import Tree
import ete3
import os

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


def GetSister(tree,target_sp):
	q_branch=tree&target_sp
	if not q_branch.get_ancestors()[0].is_root():
		for leaf in tree:
			leaf.add_features(family=leaf.name.split('|')[0])
		q_branch.add_features(family='Orobanchaceae')
		q_oro_branch=''
		#get Orobanchaceae receiver at species level
		receiver=[]
		for nd in tree.get_monophyletic(values=["Orobanchaceae"], target_attr="family"):
			if target_sp in [leaf.name for leaf in nd]:
				q_oro_branch=nd
		for leaf in q_oro_branch:
			spp=leaf.name.split('|')[-1]
			if spp.startswith('Oro') and not spp.startswith('Orobanche'):spp=spp[3:]
			try:receiver.append(id2sp[spp])
			except KeyError:receiver.append(spp)
		receiver.remove(target_sp)
		receiver.append(id2sp[target_sp.split('_')[0]])
		receiver=list(set(receiver))			
		for nd in tree.get_monophyletic(values=["Orobanchaceae"], target_attr="family"):
			if target_sp in [leaf.name for leaf in nd]:
				q_oro_branch=nd
		donor=[leaf.name for leaf in q_oro_branch.get_sisters()[0]]
		#get donor at family level
		donor_family=list(set([j.split('|')[0] for j in donor]))
		donor_genera=[j.split('|')[-1] for j in donor]
		donor_genera=list(set([j.split('_')[0] for j in donor_genera]))
		bs=tree.get_common_ancestor(donor+[target_sp]).support
		return(receiver,donor_family,donor_genera,bs)
	else:
		return('NA','NA','NA','NA')


# Function to find the largest continuous blocks with the target name value while allowing for no more than a certain number of different names (e.g., have one Asteraceae nested within 10 Orobanchaceae. This accommodates for some phylogenetic estimation error)
# Given named, ordered list
named_ordered_list = [('Nissan', 'car'), ('BMW', 'car'), ('tesla', 'car'), ('1', 'number'), ('3', 'number'), ('dodge', 'car')]

# Function to find the largest continuous block around a specific element based on its name
def find_max_block_around_element(lst, target_name, element_index, max_diff_count):
    max_block = []
    current_block = []
    diff_count = 0
    max_start = 0
    max_end = 0
    current_start = 0
    for index, item in enumerate(lst):
        name, value = item
        if index == element_index:  # When reaching the specified element index
            start_index = max(0, index - 1)  # Start searching from one element before
            end_index = min(len(lst) - 1, index + 1)  # End search at one element after
            for i in range(start_index, end_index + 1):
                current_block.append(lst[i])
                if lst[i][1] != target_name:
                    diff_count += 1
                else:
                	diff_count=0#I can have 1 Asteraceae and 1 Zingiberaceae as long as they are not next to each other 
            if diff_count <= max_diff_count and len(current_block) > len(max_block):
                max_block = current_block[:]
                max_start = start_index
                max_end = end_index
            break
    return(max_block, max_start, max_end)

# Finding the largest continuous block around a specific element based on its name
target_name = 'car'  # Replace 'car' with any target name
element_to_search = 'tesla'  # Element around which the max block is to be found
max_different_names = 1  # Maximum number of different names allowed
element_index = next(index for index, (name, _) in enumerate(named_ordered_list) if name == element_to_search)
largest_block_around_element, start_index, end_index = find_max_block_around_element(named_ordered_list, target_name, element_index, max_different_names)

# Displaying the largest continuous block around the specified element and its start and end indices
print("Largest Block around", element_to_search, ":", largest_block_around_element)
print("Start Index:", start_index)
print("End Index:", end_index)


def VGTFromTipOrder(tree,target_sp):
	tips=[node.name for node in t]
	families=[]
	orders=[]
	for j in tips:
		if j==target_sp:
			families.append('Orobanchaceae')
			orders.append('Lamiales')
		else:
			families.append(j.split('|')[0])
			if j.split('|')[0] in close_relative:orders.append('Lamiales')
			else:orders.append('NOT')
	

def HGTcalssifier(tree,target_sp):
	tips=[node.name for node in t]
	all_families=[i.split('|')[0] for i in tips]
	try:all_families.remove('NA')
	except ValueError:pass
	all_families.remove(target_sp)
	all_families=set(all_families)
	if all_families.issubset(close_relative):
		return('VGT','NA','NA','NA','NA')
	else:#use phylogeny to classify
		ancestor=t.get_midpoint_outgroup()
		t.set_outgroup(ancestor)
		q_branch=t&target_sp
		receiver,donor_family,donor_genera,bs=GetSister(t,target_sp)
		if donor_family.issubset(close_relative):return('VGT','NA','NA','NA','NA')
		else:
			#sister contains distant families
			return('HGT',';'.join(receiver),';'.join(donor_family),';'.join(donor_genera),bs)

	
targets=os.listdir('./')
targets=[i.split('.')[0] for i in targets if i.endswith('.alnmap.bed')]

for target in targets:
	blocks=open(target+'.alnmap.bed').readlines()
	out=open(target+'.hgt.sum.tsv','w')
	d=out.write('ID\tStart\tEnd\tAlignment\tClassification\tReceiver\tDonor_Family\tDonor_genera\tMethod\tBS\n')
	i=1
	for l in blocks:
		classification=''
		donor_fam=''
		donor_gen=''
		receiver=''
		method=''
		bs=''
		recs=open(target+'_HGTscanner_supporting_files/'+target+'.hgt.'+str(i)+'.fas').readlines()
		allsp=[j[1:].strip() for j in recs if j.startswith('>')]
		target_tip=[j for j in allsp if j.startswith(target)]
		allsp.remove(target_tip[0])
		oros=[j for j in allsp if j.startswith('Orobanchaceae')]
		oro_sp=[]
		genera=[j for j in allsp if not j.startswith('Orobanchaceae')]
		for j in oros:
			sp=j.split('|')[-1]
			if sp.startswith('Oro') and not sp.startswith('Orobanche'):
				sp=sp[3:]
				sp=id2sp[sp]
			oro_sp.append(sp)
		oro_sp.append(id2sp[target])
		oro_sp=list(set(oro_sp))
		families=[j.split('|')[0] for j in allsp]
		try:families.remove('NA')
		except ValueError:pass
		families=list(set(families))
		#HGT case 1: only two families (orobanchaceae and donor) are in the tree
		if len(families)<3:
			classification='High confidence HGT'
			receiver=';'.join(oro_sp)
			try:families.remove('Orobanchaceae')
			except ValueError:pass
			if len(families)<2:
				donor_fam='High confidence: '+';'.join(families)
			else:donor_fam=';'.join(families)
			genera=[j.split('|')[-1] for j in genera]
			genera=list(set([j.split('_')[0] for j in genera]))
			donor_gen=';'.join(genera)
			method='BLAST'
			bs='NA'
			d=out.write(l.strip()+'\t'+'\t'.join([classification,receiver,donor_fam,donor_gen,method,bs])+'\n')
		#HGT case 2: only one orobanchaceae is in the tree, but there are at least two other families
		elif len(oro_sp)==1:
			classification='High confidence HGT'
			method='BLAST'
			#check tree for donor
			try:
				t=Tree(target+'_HGTscanner_supporting_files/'+target+'.hgt.'+str(i)+'.aln.fas.treefile')
				reiceiver,donor_fam,donor_gen,bs=GetSister(t,oro_sp[0])
			except ete3.parser.newick.NewickError:
				receiver=target_tip[0]
				donor_fam=''
				donor_gen=''
				bs=''
			d=out.write(l.strip()+'\t'+'\t'.join([classification,';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),method,bs])+'\n')
		#use phylogeny to assess classification	
		else:
			try:
				#examine if this can be a VGT without rerooting the tree
				t=Tree(target+'_HGTscanner_supporting_files/'+target+'.hgt.'+str(i)+'.aln.fas.treefile')
				
			except ete3.parser.newick.NewickError:
				receiver=target_tip[0]
				donor_fam=''
				donor_gen=''
				bs=''
			d=out.write(l.strip()+'\t'+'\t'.join([classification,';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),method,bs])+'\n')
		i=i+1
	out.close()

