"""Plot distribution of lengths in segmented alignments."""

import matplotlib.pyplot as plt
import pandas as pd

path = r'../segment_aliscore/segment_aliscore_ungap.tsv'
df = pd.read_csv(path, sep='\t', keep_default_na=False)  # Prevents reading Asp-Ala (NA) sequences as not a number (NaN)

df_con = df[df['conserved'] == True]
df_div = df[df['conserved'] == False]

len_con = df_con.apply(lambda row: len(row['seq']), axis=1)
len_div = df_div.apply(lambda row: len(row['seq']), axis=1)

# Histogram of full range
bin_len = 10
bins_con = range(0, max(len_con) + bin_len, bin_len)
bins_div = range(0, max(len_div) + bin_len, bin_len)

plt.figure()
plt.subplots_adjust(left=0.15)  # Prevent cutoff of y-axis label
plt.hist(len_con, label='Conserved', bins=bins_con, alpha=0.5)
plt.hist(len_div, label='Diverged', bins=bins_div, alpha=0.5)
plt.title('Length Distributions of Segmented Alignment Subsequences')
plt.xlabel('Length')
plt.ylabel('Count')
plt.legend()
plt.savefig(f'len_dists_full{bin_len}.png')

# Histogram of truncated range (with smaller bins)
bin_len = 4
bins_con = range(0, max(len_con) + bin_len, bin_len)
bins_div = range(0, max(len_div) + bin_len, bin_len)

plt.figure()
plt.subplots_adjust(left=0.15)  # Prevent cutoff of y-axis label
plt.hist(len_con, label='Conserved', bins=bins_con, alpha=0.5)
plt.hist(len_div, label='Diverged', bins=bins_div, alpha=0.5)
plt.title('Length Distributions of Segmented Alignment Subsequences')
plt.xlabel('Length')
plt.ylabel('Count')
plt.legend()
plt.xlim((0, 300))
plt.savefig(f'len_dists_trun{bin_len}.png')

# Get counts
con_counts = df_con['seq'].map(lambda x: len(x)).value_counts().sort_index()
div_counts = df_div['seq'].map(lambda x: len(x)).value_counts().sort_index()

print('Conserved First Ten Counts')
print(con_counts[:10])
print()
print('Diverged First Ten Counts')
print(div_counts[:10])

"""
OUTPUT
Conserved First Ten Counts
1      92
2     226
3     222
4     206
5     247
6     297
7     373
8     270
9     168
10    282
Name: seq, dtype: int64

Diverged First Ten Counts
1      92
2     226
3     222
4     206
5     247
6     297
7     373
8     270
9     168
10    282
Name: seq, dtype: int64

DEPENDENCIES
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/segment_aliscore_ungap.tsv
"""