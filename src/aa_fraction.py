"""Calculate and plot amino acid compositions."""

import matplotlib.pyplot as plt
import pandas as pd
from numpy import arange
from sys import argv


def counts(seq):
    return {char: seq.count(char) for char in alphabet}


# Input variables
path = argv[1]  # Path to segmented sequences .tsv
key = argv[2]  # Key of column denoting subsequences class
T_name = argv[3]  # Name of True class in sentence case
F_name = argv[4]  # Name of False class in sentence case

# Constants
bar_width = 0.35
alphabet = 'DEHKRNQSTAILMVFWYCGP'

# Split into separate conserved and divergent dataframes
df = pd.read_csv(path, sep='\t', keep_default_na=False)
df_T = df[df[key] == True]
df_F = df[df[key] == False]

# Get counts as series of dict, convert series to list to dataframe, sum counts, reorder
T_counts = pd.DataFrame(df_T['seq'].apply(counts).to_list()).agg('sum')[list(alphabet)]
F_counts = pd.DataFrame(df_F['seq'].apply(counts).to_list()).agg('sum')[list(alphabet)]

# Get total sum and convert counts to fraction
T_sum = T_counts.agg('sum')
F_sum = F_counts.agg('sum')
T_fracs = T_counts / T_sum
F_fracs = F_counts / F_sum

# Plot as double bar graph
plt.figure(figsize=(8, 4))
plt.subplots_adjust(right=0.8)
index = arange(len(alphabet))
plt.bar(index, T_fracs, bar_width, align='edge', label=T_name)
plt.bar(index + bar_width, F_fracs, bar_width, align='edge', label=F_name)
plt.xlabel('Amino Acid')
plt.ylabel('Fraction')
plt.xticks(index + bar_width, T_fracs.index)
plt.title(f'Amino Acid Fractions in {T_name} and {F_name} Subsets')
plt.legend(bbox_to_anchor=(1.025, 0.5), loc="center left")
plt.savefig('aa_fraction.png')

# Print output metrics
print(f'Number of {T_name.lower()} subsequences:', len(df_T))
print(f'Number of {F_name.lower()} subsequences:', len(df_F))
print(f'Number of {T_name.lower()} amino acids:', T_sum)
print(f'Number of {F_name.lower()} amino acids:', F_sum)
