"""Plot distribution of lengths in segmented alignments."""

import matplotlib.pyplot as plt
import pandas as pd
from sys import argv

# Input variables
path = argv[1]  # Path to segmented sequences .tsv
type_name = argv[2]  # Name of column denoting segment type
T_name = argv[3]  # Name of True type in sentence case
F_name = argv[4]  # Name of False type in sentence case

# Read data and split segments
segs = pd.read_csv(path, sep='\t', keep_default_na=False)  # Prevents reading Asp-Ala (NA) sequences as not a number (NaN)
segs['seq'] = segs['seq'].map(lambda x: x.translate({ord('-'): None}))  # Remove gaps if present
T_segs = segs[segs[type_name]]
F_segs = segs[~segs[type_name]]
T_len = T_segs['seq'].map(lambda x: len(x))
F_len = F_segs['seq'].map(lambda x: len(x))

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
plt.savefig(f'len_dists_full{bin_len}.png')

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
plt.savefig(f'len_dists_trun{bin_len}.png')

# Get counts
T_counts = T_segs['seq'].map(lambda x: len(x)).value_counts().sort_index()
F_counts = F_segs['seq'].map(lambda x: len(x)).value_counts().sort_index()

print(f'{T_name} First Ten Counts')
print(T_counts[:10])
print()
print(f'{F_name} First Ten Counts')
print(F_counts[:10])
