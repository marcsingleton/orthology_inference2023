"""Plot distribution of lengths in segmented alignments."""

import matplotlib.pyplot as plt
import pandas as pd
from sys import argv

# Input variables
path = argv[1]  # Path to segmented sequences .tsv
key = argv[2]  # Key of column denoting subsequences class
T_name = argv[3]  # Name of True class in sentence case
F_name = argv[4]  # Name of False class in sentence case

# Read data and split subsequences
df = pd.read_csv(path, sep='\t', keep_default_na=False)  # Prevents reading Asp-Ala (NA) sequences as not a number (NaN)
df['seq'] = df['seq'].map(lambda x: x.translate({ord('-'): None}))  # Remove gaps if present
df_T = df[df[key] == True]
df_F = df[df[key] == False]
len_T = df_T['seq'].map(lambda x: len(x))
len_F = df_F['seq'].map(lambda x: len(x))

# Histogram of full range
bin_len = 10
bins_T = range(0, max(len_T) + bin_len, bin_len)
bins_F = range(0, max(len_F) + bin_len, bin_len)

plt.figure()
plt.subplots_adjust(left=0.15)  # Prevent cutoff of y-axis label
plt.hist(len_T, label=T_name, bins=bins_T, alpha=0.5)
plt.hist(len_F, label=F_name, bins=bins_F, alpha=0.5)
plt.title('Length Distributions of Segmented Alignment Subsequences')
plt.xlabel('Length')
plt.ylabel('Count')
plt.legend()
plt.savefig(f'len_dists_full{bin_len}.png')

# Histogram of truncated range (with smaller bins)
bin_len = 4
bins_T = range(0, max(len_T) + bin_len, bin_len)
bins_F = range(0, max(len_F) + bin_len, bin_len)

plt.figure()
plt.subplots_adjust(left=0.15)  # Prevent cutoff of y-axis label
plt.hist(len_T, label=T_name, bins=bins_T, alpha=0.5)
plt.hist(len_F, label=F_name, bins=bins_F, alpha=0.5)
plt.title('Length Distributions of Segmented Alignment Subsequences')
plt.xlabel('Length')
plt.ylabel('Count')
plt.legend()
plt.xlim((0, 300))
plt.savefig(f'len_dists_trun{bin_len}.png')

# Get counts
T_counts = df_T['seq'].map(lambda x: len(x)).value_counts().sort_index()
F_counts = df_F['seq'].map(lambda x: len(x)).value_counts().sort_index()

print(f'{T_name} First Ten Counts')
print(T_counts[:10])
print()
print(f'{F_name} First Ten Counts')
print(F_counts[:10])
