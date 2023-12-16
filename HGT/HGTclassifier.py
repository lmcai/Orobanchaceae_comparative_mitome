from ete3 import Tree
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


def HGTcalssifier(treefile,target_sp):
	t=Tree(treefile)
	tips=[node.name for node in t]
	all_families=[i.split('|')[0] for i in tips]
	try:all_families.remove('NA')
	except ValueError:pass
	all_families.remove(target_sp)
	all_families=set(all_families)
	if all_families.issubset(close_relative):return('VGT')
	else:#use phylogeny to classify
		ancestor=t.get_midpoint_outgroup()
		t.set_outgroup(ancestor)
		q_branch=t&target_sp
		if not q_branch.get_ancestors()[0].is_root():
			for leaf in t:
				leaf.add_features(family=leaf.name.split('|')[0])
			q_branch.add_features(family='Orobanchaceae')
			q_oro_branch=''
			#get Orobanchaceae receiver at species level
			receiver=[]
			for nd in t.get_monophyletic(values=["Orobanchaceae"], target_attr="family"):
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
			donor=[leaf.name for leaf in q_oro_branch.get_sisters()[0]]
			#get donor at family level
			donor_family=set([j.split('|')[0] for j in donor])
			if donor_family.issubset(close_relative):return('VGT')
			else:
				#sister contains distant families
				return(receiver,donor_family,donor_genera)
		else:
			return('')
	
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
			classification='HGT'
			receiver=';'.join(oro_sp)
			try:families.remove('Orobanchaceae')
			except ValueError:pass
			donor_fam=';'.join(families)
			genera=[j.split('|')[-1] for j in genera]
			genera=list(set([j.split('_')[0] for j in genera]))
			donor_gen=';'.join(genera)
			method='BLAST'
			bs='na'
			d=out.write(l.strip()+'\t'+'\t'.join([classification,receiver,donor_fam,donor_gen,method,bs])+'\n')
		#HGT case 2: only one orobanchaceae is in the tree
		elif len(oro_sp)==1:
			classification='HGT'
			receiver=oro_sp[0]
			method='BLAST'
			#check tree for donor
			donor_fam=''
			donor_gen=''
			bs=''
			d=out.write(l.strip()+'\t'+'\t'.join([classification,receiver,donor_fam,donor_gen,method,bs])+'\n')
		else:
			#need phylogeny for classification	
			pass
		i=i+1
	
	out.close()

