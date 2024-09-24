import json
from ete3 import Tree

with open('all_mt_no_ribosome.conc.fasta.RELAX.json', 'r') as file:
	data = json.load(file)

k_values = {name: info["k (general descriptive)"] for name, info in data["branch attributes"]["0"].items()}

phylo=data["input"]["trees"]["0"]
phylo=Tree(phylo+';',format=1)

#annotate k values as branch lengths
for node in phylo.traverse("postorder"):
	if node.name in k_values.keys():
		node.dist=k_values[node.name]

	
phylo.write(format=1)