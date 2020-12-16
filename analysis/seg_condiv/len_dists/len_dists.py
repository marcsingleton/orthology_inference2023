"""Plot distribution of lengths in segmented alignments."""

import matplotlib.pyplot as plt
import os
import pandas as pd

# Input variables
path = '../segment_aliscore/out/segment_aliscore.tsv'  # Path to segmented sequences .tsv
type_name = 'conserved'  # Name of column denoting segment type
T_name = 'Conserved'  # Name of True type in sentence case
F_name = 'Diverged'  # Name of False type in sentence case

# Read data and split segments
segs = pd.read_csv(path, sep='\t', keep_default_na=False)  # Prevents reading Asp-Ala (NA) sequences as not a number (NaN)
segs['seq'] = segs['seq'].map(lambda x: x.translate({ord('-'): None}))  # Remove gaps if present
T_segs = segs[segs[type_name]]
F_segs = segs[~segs[type_name]]
T_len = T_segs['seq'].map(lambda x: len(x))
F_len = F_segs['seq'].map(lambda x: len(x))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Histogram of full range
bin_len = 10
T_bins = range(0, max(T_len) + bin_len, bin_len)
F_bins = range(0, max(F_len) + bin_len, bin_len)

plt.figure()
plt.subplots_adjust(left=0.15)  # Prevent cutoff of y-axis label
plt.hist(T_len, label=T_name, bins=T_bins, alpha=0.5)
plt.hist(F_len, label=F_name, bins=F_bins, alpha=0.5)
plt.title('Length Distributions of Segmented Alignment Subsequences')
plt.xlabel('Length')
plt.ylabel('Count')
plt.legend()
plt.savefig(f'out/len_dists_full{bin_len}.png')

# Histogram of truncated range (with smaller bins)
bin_len = 4
T_bins = range(0, max(T_len) + bin_len, bin_len)
F_bins = range(0, max(F_len) + bin_len, bin_len)

plt.figure()
plt.subplots_adjust(left=0.15)  # Prevent cutoff of y-axis label
plt.hist(T_len, label=T_name, bins=T_bins, alpha=0.5)
plt.hist(F_len, label=F_name, bins=F_bins, alpha=0.5)
plt.title('Length Distributions of Segmented Alignment Subsequences')
plt.xlabel('Length')
plt.ylabel('Count')
plt.legend()
plt.xlim((0, 300))
plt.savefig(f'out/len_dists_trun{bin_len}.png')

# Get counts
T_counts = T_segs['seq'].map(lambda x: len(x)).value_counts().sort_index()
F_counts = F_segs['seq'].map(lambda x: len(x)).value_counts().sort_index()

print(f'{T_name} First Ten Counts')
print(T_counts[:10])
print()
print(f'{F_name} First Ten Counts')
print(F_counts[:10])

"""
OUTPUT
Conserved First Ten Counts
0    9900
1    3303
2    5834
3    7384
4    7714
5    7937
6    7674
7    6956
8    6492
9    5629
Name: seq, dtype: int64

Diverged First Ten Counts
0    70061
1    33157
2    29108
3    24995
4    20988
5    17086
6    14100
7    11946
8     9917
9     8073
Name: seq, dtype: int64

DEPENDENCIES
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/out/segment_aliscore.tsv
"""