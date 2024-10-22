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
"Lcl":"Lathraea_clandestina",
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
"Lcl":"Lathraea_clandestina",
"Lsq":"Lathraea_squamaria"
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
	#extend on the left side
	for i in range(current_start,0,-1):
		tip, taxonomy = lst[i]
		if taxonomy == taxonomy_name:  # When same as the target taxon rank
			current_start = max(0, current_start - 1)  # Start searching from one element before
			diff_count=0
			first_ingroup=i
			ingroup_sp.append(tip)
		else:
			if diff_count< max_diff_count:current_start = max(0, current_start - 1)
			else:break
			diff_count += 1
	#extend on the right side
	for i in range(current_end,len(lst)):
		tip, taxonomy = lst[i]
		if taxonomy == taxonomy_name:  # When same as the target taxon rank
			current_end = min(len(lst), current_end+1)  # Start searching from one element before
			diff_count=0
			last_ingroup=i
			ingroup_sp.append(tip)
		else:
			if diff_count< max_diff_count:current_end= min(len(lst), current_end+1)
			else:break
			diff_count += 1
	return(first_ingroup, last_ingroup, ingroup_sp)

#take a tree and return true or false for VGT based on tip orders, this is a generous way to classify VGT and better accommodates topology error in a phylo tree
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
	max_different_names = 1  # Maximum number of different names allowed
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


def GetSister(tree,target_sp):
	q_branch=tree&target_sp
	if q_branch.get_ancestors()[0].is_root():
		ancestor=tree.get_midpoint_outgroup()
		tree.set_outgroup(ancestor)
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
			if spp==target_sp:
				receiver.append(id2sp[target_sp.split('_')[0]])
			elif spp.startswith('Oro') and not spp.startswith('Orobanche'):
				spp=spp[3:]
				try:receiver.append(id2sp[spp])
				except KeyError:receiver.append(spp)
			else:receiver.append(spp)
		receiver=list(set(receiver))			
		for nd in tree.get_monophyletic(values=["Orobanchaceae"], target_attr="family"):
			if target_sp in [leaf.name for leaf in nd]:
				q_oro_branch=nd
		donor=[leaf.name for leaf in q_oro_branch.get_sisters()[0]]
		#get donor at family level
		donor_family=list(set([j.split('|')[0] for j in donor]))
		donor_genera_temp=[j.split('|')[-1] for j in donor]
		donor_genera_temp=list(set([j.split('_')[0] for j in donor_genera_temp]))
		donor_genera=[]
		for k in donor_genera_temp:
			if k.startswith('Oro') and not k.startswith('Orobanche'):
				try:donor_genera.append(id2sp[k[3:]].split('_')[0])
				except KeyError:donor_genera.append(k)
			else:donor_genera.append(k)
		donor_genera=list(set(donor_genera))
		bs=tree.get_common_ancestor(donor+[target_sp]).support
		return(receiver,donor_family,donor_genera,str(bs))
	else:
		return('NA','NA','NA','NA')

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
		donor_family=set(donor_family)
		if donor_family.issubset(close_relative):return('VGT','NA','NA','NA','NA')
		else:
			#sister contains distant families
			return('HGT',';'.join(receiver),';'.join(donor_family),';'.join(donor_genera),bs)


def main(target):
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
		oros.append(target_tip[0])
		oro_sp=ID2sp(oros,target_tip[0])
		genera=[j for j in allsp if not j.startswith('Orobanchaceae')]
		families=[j.split('|')[0] for j in allsp]
		families=list(set(families))
		try:families.remove('Orobanchaceae')
		except ValueError:pass
		if not families==['NA']:
			try:families.remove('NA')
			except ValueError:pass
		#HGT case 1: only two families (orobanchaceae and donor) are in the tree
		if len(families)<2:
			if families[0] in close_relative:
				#A VGT case since the only other family is a close relative
				d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','BLAST','NA'])+'\n')
			else:
				classification='High confidence HGT'
				receiver=';'.join(oro_sp)
				donor_fam='High confidence: '+';'.join(families)
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
				receiver,donor_fam,donor_gen,bs=GetSister(t,target_tip[0])
				d=out.write(l.strip()+'\t'+'\t'.join([classification,';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),method,bs])+'\n')
			except ete3.parser.newick.NewickError:
				receiver=target_tip[0]
				donor_fam=''
				donor_gen=''
				bs=''
				if len(allsp)<4:
					#too few tips to run a tree
					donor_fam=families
					genera=[j.split('|')[-1] for j in genera]
					genera=list(set([j.split('_')[0] for j in genera]))
					donor_gen=genera
					bs='NA'
				d=out.write(l.strip()+'\t'+'\t'.join([classification,receiver,';'.join(donor_fam),';'.join(donor_gen),method,bs])+'\n')
		#use phylogeny to assess classification	
		else:
			try:
				method='Phylogeny'
				#examine if this can be a VGT without rerooting the tree
				t=Tree(target+'_HGTscanner_supporting_files/'+target+'.hgt.'+str(i)+'.aln.fas.treefile')
				#VGT based on tip order
				if VGTFromTipOrder(t,target_tip[0]):
					d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','Phylogeny','NA'])+'\n')
				else:
					#use tree topology to make decisions
					receiver,donor_fam,donor_gen,bs=GetSister(t,target_tip[0])
					if len(donor_fam)==1:
						if not donor_fam[0] in close_relative and float(bs)>85:
							classification='High confidence HGT'
							d=out.write(l.strip()+'\t'+'\t'.join([classification,';'.join(receiver),'High confidence: '+';'.join(donor_fam),';'.join(donor_gen),method,bs])+'\n')
						else:
							if donor_fam[0] in close_relative:
								#this family is within lamiales
								d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','Phylogeny','NA'])+'\n')
							else:
								#one family that's not close relative, but with low support
								classification='Putative HGT'
								d=out.write(l.strip()+'\t'+'\t'.join([classification,';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),method,bs])+'\n')
					else:
						#multiple donor families
						classification='Inconclusive'
						d=out.write(l.strip()+'\t'+'\t'.join([classification,';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),method,bs])+'\n')
			except ete3.parser.newick.NewickError:
				receiver=target_tip[0]
				donor_fam=''
				donor_gen=''
				bs=''
				d=out.write(l.strip()+'\t'+'\t'.join([classification,receiver,';'.join(donor_fam),';'.join(donor_gen),method,bs])+'\n')
		i=i+1
	out.close()

main(sys.argv[1])