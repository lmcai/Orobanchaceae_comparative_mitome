from ete3 import Tree

close_relative={'Orobanchaceae','Lamiaceae','Bignoniaceae','Phrymaceae','Acanthaceae','Oleaceae','Scrophulariaceae','Verbenaceae','Plantaginaceae','Gesneriaceae','Lentibulariaceae'}
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
	all_families.remove('NA')
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
			#get Orobanchaceae receiver at genus level
			q_oro_branch=''
			receiver=[]
			for nd in t.get_monophyletic(values=["Orobanchaceae"], target_attr="family"):
				if target_sp in [leaf.name for leaf in nd]:
					q_oro_branch=nd
			for leaf in q_oro_branch:
				try:receiver.append(id2sp[leaf.name.split('|')[-1]])
				except KeyError:
					spp=leaf.name.split('|')[-1]
					spp=spp[3:]
					receiver.append(spp)
			receiver=list(set([j.split('_')[0] for j in receiver]))
			sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
		else:
			print(treefile,'bad rooting')
	

Evaluate the source of the region and output summary file
out=open(sp+'.hgt.sum.tsv','w')
out.write('ID\tTarget_scaffold\tStart\tEnd\tPhylo_source\tBlast_hit_ID\n')
out.close()
out=open(sp+'.hgt.sum.tsv','a')

#for i in range(1,order):
#	q=loci[i-1].split()[0]
#	outgroup=[]
#	try:
#		t=Tree(sp+'.gt.'+str(i)+'.aln.fas.treefile')
#		ancestor=t.get_midpoint_outgroup()
#		t.set_outgroup(ancestor)
#		q_branch=t&q
#		if not q_branch.get_ancestors()[0].is_root():
#			sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
#		else:
#			t=Tree(sp+'.gt.'+str(i)+'.aln.fas.treefile')
#			outgroup=[leaf.name for leaf in t if leaf.name.startswith('Sorghum')]
#			if len(outgroup)>0:
#				t.set_outgroup(outgroup[0])
#				q_branch=t&q
#				sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
#			else:
#				outgroup=[leaf.name for leaf in t if leaf.name.startswith('Rehm')]
#				if outgroup:
#					t.set_outgroup(t&outgroup[0])
#					q_branch=t&q
#					sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
#				else:
					#no sorghum no rehmannia
#					sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
#		out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
#	except ete3.parser.newick.NewickError:
#		sisters=open(sp+'.temp.'+str(i)+".fas").readlines()
#		sisters=[i[1:].strip() for i in sisters if (i.startswith('>')) and (not i[1:].strip()==q)]
#		out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+'ALL HITS: '+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
		
#out.close()
