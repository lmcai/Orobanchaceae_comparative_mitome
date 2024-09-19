import matplotlib.pyplot as plt
import sys

pos=open(sys.argv[1]).readlines()

HGT_dist=[]
cur_loci=pos[1].split()[0]
cur_end='NA'
for l in pos[1:]:
    if l.split()[0] == cur_loci:
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