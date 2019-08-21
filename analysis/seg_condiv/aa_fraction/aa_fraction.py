"""Calculate and plot amino acid compositions."""

import matplotlib.pyplot as plt
import pandas as pd
from numpy import arange


def counts(seq):
    return {char: seq.count(char) for char in alphabet}


path = '../segment_aliscore/segment_aliscore_ungap.tsv'
bar_width = 0.35
alphabet = 'DEHKRNQSTAILMVFWYCGP'


# Split into separate conserved and divergent dataframes
df = pd.read_csv(path, sep='\t', keep_default_na=False)
df_con = df[df['conserved'] == True]
df_div = df[df['conserved'] == False]

# Get counts as series of dict, convert series to list to dataframe, sum counts, reorder
con_counts = pd.DataFrame(df_con['seq'].apply(counts).to_list()).agg('sum')[list(alphabet)]
div_counts = pd.DataFrame(df_div['seq'].apply(counts).to_list()).agg('sum')[list(alphabet)]

# Get total sum and convert counts to fraction
con_sum = con_counts.agg('sum')
div_sum = div_counts.agg('sum')
con_fracs = con_counts / con_sum
div_fracs = div_counts / div_sum

# Plot as double bar graph
plt.figure(figsize=(8, 4))
plt.subplots_adjust(right=0.8)
index = arange(len(alphabet))
plt.bar(index, con_fracs, bar_width, align='edge', label='Conserved')
plt.bar(index + bar_width, div_fracs, bar_width, align='edge', label='Diverged')
plt.xlabel('Amino Acid')
plt.ylabel('Fraction')
plt.xticks(index + bar_width, con_fracs.index)
plt.title('Amino Acid Fractions in Conserved and Divergent Subsets')
plt.legend(bbox_to_anchor=(1.025, 0.5), loc="center left")
plt.savefig('aa_fraction.png')

# Print output metrics
print('Number of conserved subsequences:', len(df_con))
print('Number of diverged subsequences:', len(df_div))
print('Number of conserved amino acids:', con_sum)
print('Number of diverged amino acids:', div_sum)

"""
OUTPUT
Number of conserved subsequences: 322940
Number of diverged subsequences: 305170
Number of conserved amino acids: 25887032
Number of diverged amino acids: 2197778

DEPENDENCIES
../segment_aliscore/segment_aliscore_ungap.py
    ../segment_aliscore/segment_aliscore_ungap.tsv
"""