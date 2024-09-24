import json
import ete3

file=open('all_mt_no_ribosome.hyphy2.5.8.2rate.RELAX.json', 'r')
data = json.load(file)

k_values = {}
#omega={}
general_des_omega={}
# Navigate to the "branch attributes" section
branch_attributes = data['branch attributes']['0']

# Iterate over each species and extract the k value
for species, attributes in branch_attributes.items():
    k_values[species] = attributes.get("k (general descriptive)")
    #omega[species] = attributes.get("MG94xREV with separate rates for branch sets")
    general_des_omega[species] = attributes.get("General descriptive")


##########################
#output to spreadsheet
out=open('mt_no_ribosome.hyphy.param.tsv','w')
out.write('Species/node\tk\tomega\n')
for k in k_values.keys():
	d=out.write(k+'\t'+str(k_values[k])+'\t'+str(general_des_omega[k])+'\n')
out.close()

#############################
#tree annotation
t=data['input']['trees']['0']
t=Tree(t,format=1)

for node in t.traverse():
	if node.is_leaf():
		pass
		#if node.name in k_values.keys():node.add_features(k=k_values[node.name])
	else:
		if node.name in k_values.keys():node.add_features(k=k_values[node.name])
    

t.write(format=7,features=["k"])
