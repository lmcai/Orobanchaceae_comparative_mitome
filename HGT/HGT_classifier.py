from ete3 import Tree

close_relative={'Orobanchaceae','Lamiaceae','Bignoniaceae','Phrymaceae','Acanthaceae','Oleaceae','Scrophulariaceae','Verbenaceae','Plantaginaceae','Gesneriaceae','Lentibulariaceae','Apocynaceae'}
id2sp={'Ain':"Aeginetia_indica",
"Ase":"Alectra_sessiflora",
"Afa":"Aphyllon_faciculatum",
"Apu":"Aphyllon_purpureum",
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


def HGT_calssifier(treefile,target_sp):
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
			#get Orobanchaceae receiver at genus level
			receiver=[]
			for nd in t.get_monophyletic(values=["Orobanchaceae"], target_attr="family"):
				if target_sp in [leaf.name for leaf in nd]:
					q_oro_branch=nd
			for leaf in q_oro_branch:
				spp=leaf.name.split('|')[-1]
				if spp.startswith('Oro'):spp=spp[3:]
				try:receiver.append(id2sp[spp])
				except KeyError:receiver.append(spp)
			receiver.remove(target_sp)
			receiver.append(id2sp[target_sp.split('_')[0]])
			receiver_genus=list(set([j.split('_')[0] for j in receiver]))
			donor=[leaf.name for leaf in q_oro_branch.get_sisters()[0]]
			#get donor at family level
			donor_family=set([j.split('|')[0] for j in donor])
			if donor_family.issubset(close_relative):return('VGT')
			else:
				#sister contains distant families
				print(donor_family)
		else:
			print(treefile,'bad rooting')
	

#Evaluate the source of the region and output summary file
out=open(sp+'.hgt.sum.tsv','w')
out.write('ID\tTarget_scaffold\tStart\tEnd\tPhylo_source\tBlast_hit_ID\n')
out.close()
out=open(sp+'.hgt.sum.tsv','a')
