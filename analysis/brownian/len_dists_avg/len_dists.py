"""Plot distribution of lengths in segmented alignments."""

import matplotlib.pyplot as plt
import os
import pandas as pd

# Input variables
path = '../segment_avg/out/segment_avg.tsv'  # Path to segmented sequences .tsv
type_name = 'ordered'  # Name of column denoting segment type
T_name = 'Ordered'  # Name of True type in sentence case
F_name = 'Disordered'  # Name of False type in sentence case

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
Ordered First Ten Counts
0    20760
1     4118
2     5636
3     6614
4     6671
5     6452
6     6354
7     5346
8     4991
9     4420
Name: seq, dtype: int64

Disordered First Ten Counts
0    20944
1     5414
2     8450
3     8917
4     9558
5     8761
6     8357
7     7315
8     6521
9     5705
Name: seq, dtype: int64

DEPENDENCIES
../segment_avg/segment_avg.py
    ../segment_avg/out/segment_avg.tsv
"""