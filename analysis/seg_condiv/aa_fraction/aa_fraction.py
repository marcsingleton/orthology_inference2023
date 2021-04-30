"""Calculate and plot amino acid compositions."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from numpy import arange


def counts(seq):
    return {char: seq.count(char) for char in alphabet}


# Input variables
path = '../segment_aliscore/out/segment_aliscore.tsv'  # Path to segmented sequences .tsv
type_name = 'conserved'  # Name of column denoting segment type
T_name = 'Conserved'  # Name of True type in sentence case
F_name = 'Diverged'  # Name of False type in sentence case

# Constants
bar_width = 0.35
alphabet = 'DEHKRNQSTAILMVFWYCGP'

# Read data and split segments
segs = pd.read_csv(path, sep='\t', keep_default_na=False)
T_segs = segs[segs[type_name]]
F_segs = segs[~segs[type_name]]

# Get counts as series of dict, convert series to list to dataframe, sum counts, reorder
T_counts = pd.DataFrame(T_segs['seq'].apply(counts).to_list()).agg('sum')[list(alphabet)]
F_counts = pd.DataFrame(F_segs['seq'].apply(counts).to_list()).agg('sum')[list(alphabet)]

# Get total sum and convert counts to fraction
T_sum = T_counts.agg('sum')
F_sum = F_counts.agg('sum')
T_fracs = T_counts / T_sum
F_fracs = F_counts / F_sum

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

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
plt.legend(bbox_to_anchor=(1.025, 0.5), loc='center left')
plt.savefig('out/aa_fraction.png')

# Print output metrics
print(f'Number of {T_name.lower()} subsequences:', len(T_segs))
print(f'Number of {F_name.lower()} subsequences:', len(F_segs))
print(f'Number of {T_name.lower()} amino acids:', T_sum)
print(f'Number of {F_name.lower()} amino acids:', F_sum)

"""
OUTPUT
Number of conserved subsequences: 314480
Number of diverged subsequences: 296370
Number of conserved amino acids: 26124728
Number of diverged amino acids: 1960082

DEPENDENCIES
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/out/segment_aliscore.tsv
"""