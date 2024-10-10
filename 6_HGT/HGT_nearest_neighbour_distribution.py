import matplotlib.pyplot as plt
import sys
import numpy as np

#usage:
#python HGT_nearest_neighbour_distribution.py Afa.hgt.sum.tsv
pos=open(sys.argv[1]).readlines()

HGT_dist=[]
cur_loci=pos[1].split()[0]
cur_end='NA'
for l in pos[1:]:
	if l.split()[0] == cur_loci:
		#if l.split('\t')[4]=='High confidence HGT':
		if 'HGT' in l.split('\t')[4]:
			if cur_end == 'NA':
				cur_end = int(l.split()[2])
			else:
				# Append the difference between the first part and cur_end to HGT_dist
				if (int(l.split()[1]) - cur_end)>0:HGT_dist.append(int(l.split()[1]) - cur_end)
				# Update cur_end to the second part of the current element
				cur_end = int(l.split()[2])    
		else:
			# If 'HGT' is not in the fifth part, do nothing
			pass
	else:
		# Update cur_loci to the first part of the current element
		cur_loci = l.split()[0]
		cur_end='NA'
		#if l.split('\t')[4]=='High confidence HGT':
		if 'HGT' in l.split('\t')[4]:
			cur_end = int(l.split()[2])

# Plot histogram
plt.hist(HGT_dist, bins=50, color='skyblue', edgecolor='black')

# Add title
#plt.title('Histogram of Integers')

# Add length and width
plt.xlabel('Distance to nearest neighbour (bp)')
plt.ylabel('Frequency')

# Adjust layout to prevent clipping of labels
plt.tight_layout()

# Save to PDF
plt.savefig(sys.argv[1].split('.')[0]+'.hgt_neighbour_dist.pdf')

# Show plot (optional)
#plt.show()

############################
#test for random distribution
size_file=open('sp_mito_hgt_size.csv').readlines()
number_hgt=len(HGT_dist)
sp=sys.argv[1].split('.')[0]
for l in size_file:
	if l.split(',')[0]==sp:genome_length = int(l.split(',')[1])-int(l.split(',')[2].strip())

def nearest_neighbor_distances(positions):
    distances = np.diff(positions)
    return distances

def simulate_null_distribution(num_genes, genome_length, num_simulations=1000):
    null_means = []
    null_variances = []
    for _ in range(num_simulations):
        random_positions = np.sort(np.random.randint(0, genome_length, num_genes))
        distances = nearest_neighbor_distances(random_positions)
        mean_dist, var_dist = calculate_mean_variance(distances)
        null_means.append(mean_dist)
        null_variances.append(var_dist)
    return np.array(null_means), np.array(null_variances)

def calculate_mean_variance(distances):
    mean_dist = np.mean(distances)
    var_dist = np.var(distances)
    return mean_dist, var_dist

observed_mean, observed_variance = calculate_mean_variance(HGT_dist)

#Simulate null distribution
null_means, null_variances = simulate_null_distribution(number_hgt, genome_length)

#Calculate p-values
p_value_mean = np.mean(null_means <= observed_mean)
p_value_variance = np.mean(null_variances <= observed_variance)

#print(f"Observed mean: {observed_mean}, p-value (mean): {p_value_mean}")
#print(f"Observed variance: {observed_variance}, p-value (variance): {p_value_variance}")

# Optional: Combine p-values using Fisher's method
from scipy.stats import chi2
combined_chi2 = -2 * (np.log(p_value_mean) + np.log(p_value_variance))
combined_p_value = chi2.sf(combined_chi2, 4)  # Degrees of freedom = 4
print(f"{sp} Combined p-value: {combined_p_value}")
