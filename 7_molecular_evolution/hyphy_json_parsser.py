import json
import ete3

file=open('all_mt_no_ribosome.hyphy2.5.8.2rate.RELAX.json', 'r')
data = json.load(file)

k_values = {}

# Navigate to the "branch attributes" section
branch_attributes = data['branch attributes']['0']

# Iterate over each species and extract the k value
for species, attributes in branch_attributes.items():
    k_value = attributes.get("k (general descriptive)")
    k_values[species] = k_value

t=data['input']['trees']['0']
t=Tree(t,format=1)

for node in t.traverse():
	if node.is_leaf():
		pass
		#if node.name in k_values.keys():node.add_features(k=k_values[node.name])
	else:
		if node.name in k_values.keys():node.add_features(k=k_values[node.name])
    

t.write(format=7,features=["k"])
