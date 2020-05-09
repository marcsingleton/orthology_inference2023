"""Plot various statistics of OGs."""

import matplotlib.pyplot as plt
import os
import pandas as pd

# Load rOGs
df = pd.read_csv('../reduce_gclusters/out/rclusters.tsv', sep='\t')
groups = df.groupby('rOGid')
num_seqs = df['ppid'].nunique()
num_rOGs = df['rOGid'].nunique()

# Make output directory
if not os.path.exists(f'out/'):
    os.mkdir(f'out/')

# Plot stats
# Representation of species in rOGs
spid_rOGs = {}
for _, group in groups:
    for spid in group['spid'].unique():
        spid_rOGs[spid] = spid_rOGs.get(spid, 0) + 1

labels, h_rOG = zip(*sorted(spid_rOGs.items(), key=lambda i: i[0]))
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_rOG, tick_label=labels, align='center')
ax1.set_xlabel('Species')
ax1.set_ylabel('Count')
ax1.set_title('Representation of Species in rOGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_rOGs, mx / num_rOGs)
ax2.set_ylabel('Fraction')

fig.tight_layout()
fig.savefig('out/species_rOG_rep.png')
plt.close()

# Distribution of species in sequences
spid_seqs = df['spid'].value_counts().sort_index()
labels, h_seq = zip(*spid_seqs.items())
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_seq, tick_label=labels)
ax1.set_xlabel('Species')
ax1.set_ylabel('Count')
ax1.set_title('Distribution of Species in Sequences')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_seqs, mx / num_seqs)
ax2.set_ylabel('Fraction')

fig.tight_layout()
fig.savefig('out/species_seq_dist.png')
plt.close()

# Correlation of rOG and sequence counts
fig, ax = plt.subplots()
ax.scatter(h_rOG, h_seq)
ax.set_xlabel('Species Count in rOGs')
ax.set_ylabel('Species Count in Sequences')
ax.set_title('rOG and Sequence Count Correlation')

fig.savefig('out/species_rOGseq_corr.png')
plt.close()

# Distribution of number of species
dist_species = groups['spid'].nunique().value_counts()
spec, spec_count = zip(*dist_species.items())
fig, ax1 = plt.subplots()
ax1.bar(spec, spec_count)
ax1.set_title('Distribution of Number of Species in rOGs')
ax1.set_xlabel('Number of Species')
ax1.set_ylabel('Count')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_rOGs, mx / num_rOGs)
ax2.set_ylabel('Fraction')

fig.tight_layout()
fig.savefig('out/num_species_dist.png')
plt.close()

# Distribution of number of sequences
dist_seq = groups.size().value_counts()
seq, seq_count = zip(*dist_seq.items())
fig, ax1 = plt.subplots()
ax1.bar(seq, seq_count, width=1, align='center')
ax1.set_title('Distribution of Number of Sequences in rOGs')
ax1.set_xlabel('Number of Sequences')
ax1.set_ylabel('Count')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_rOGs, mx / num_rOGs)
ax2.set_ylabel('Fraction')

fig.tight_layout()
fig.savefig('out/num_sequence_dist.png')
plt.close()

# Distribution of number of duplicates
dist_dup = (groups.size() - groups['spid'].nunique()).value_counts()
seq, seq_count = zip(*dist_dup.drop(0, errors='ignore').items())  # Drop 0; ignore error if 0 does not exist
fig, ax1 = plt.subplots()
ax1.bar(seq, seq_count, width=1, align='center')
ax1.set_title('Distribution of Number of Duplicates in rOGs')
ax1.set_xlabel('Number of Sequences')
ax1.set_ylabel('Count')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_rOGs, mx / num_rOGs)
ax2.set_ylabel('Fraction')

fig.tight_layout()
fig.savefig('out/num_duplicate_dist.png')
plt.close()

# Print counts
print('number of alignments:', num_rOGs)
print()
print('number of alignments with 10 species:', dist_species[10])
print('fraction of alignments with 10 species:', dist_species[10] / num_rOGs)
print()
print('number of alignments with 10 sequences:', dist_seq[10])
print('fraction of alignments with 10 sequences:', dist_seq[10] / num_rOGs)
print()
print('number of alignments with duplicates:', num_rOGs - dist_dup[0])
print('fraction of alignments with duplicates', (num_rOGs - dist_dup[0]) / num_rOGs)

"""
OUTPUT
number of alignments: 12885

number of alignments with 10 species: 10182
fraction of alignments with 10 species: 0.790221187427241

number of alignments with 10 sequences: 9074
fraction of alignments with 10 sequences: 0.7042297244858362

number of alignments with duplicates: 1368
fraction of alignments with duplicates 0.10616996507566938

NOTES
These plots are largely based off those in analysis/EggNOGv5_validation/ali_stats/ali_stats.py

DEPENDENCIES
../reduce_gclusters/reduce_gclusters.py
    ../reduce_gclusters/out/rclusters.tsv
"""