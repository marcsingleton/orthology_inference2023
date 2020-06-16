"""Plot various statistics of OGs."""

import matplotlib.pyplot as plt
import os
import pandas as pd

# Load gn metadata
gnid2spid = {}
with open('../ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        _, gnid, spid = line.split()
        gnid2spid[gnid] = spid

# Load OGs
rows = []
with open('../cluster_xgraph/out/gclusters.tsv') as file:
    for i, line in enumerate(file):
        OGid = hex(i)[2:].zfill(4)
        gnids = set([node for edge in line.split() for node in edge.split(',')])
        for gnid in gnids:
            rows.append({'OGid': OGid, 'gnid': gnid, 'spid': gnid2spid[gnid]})
df = pd.DataFrame(rows)

groups = df.groupby('OGid')
num_OGs = df['OGid'].nunique()
num_ugns = df['gnid'].nunique()
num_u10 = len(groups.filter(lambda x: len(x) == 10 and x['spid'].nunique() == 10)) // 10

# Make output directory
if not os.path.exists(f'out/'):
    os.mkdir(f'out/')

# Plot stats
# Number of associated OGs for species
spid_OGs = {}
for _, group in groups:
    for spid in group['spid'].unique():
        spid_OGs[spid] = spid_OGs.get(spid, 0) + 1

labels, h_OG = zip(*sorted(spid_OGs.items(), key=lambda i: i[0]))
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_OG, tick_label=labels, align='center')
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of Associated OGs')
ax1.set_title('Number of Associated OGs for each Species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_OGs, mx / num_OGs)
ax2.set_ylabel('Fraction of Total OGs')

fig.tight_layout()
fig.savefig('out/OGnum-species.png')
plt.close()

# Distribution of genes across species
spid_gns = df['spid'].value_counts().sort_index()
labels, h_gn = zip(*spid_gns.items())
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_gn, tick_label=labels)
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of Genes')
ax1.set_title('Distribution of Genes across Species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / len(df), mx / len(df))
ax2.set_ylabel('Fraction of Total Genes')

fig.tight_layout()
fig.savefig('out/gnnum-species_dist.png')
plt.close()

# Distribution of unique genes across species
spid_ugns = df.loc[:, ['gnid', 'spid']].drop_duplicates()['spid'].value_counts().sort_index()
labels, h_ugn = zip(*spid_ugns.items())
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_ugn, tick_label=labels)
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of Unique Genes')
ax1.set_title('Distribution of Unique Genes across Species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_ugns, mx / num_ugns)
ax2.set_ylabel('Fraction of Total Unique Genes')

fig.tight_layout()
fig.savefig('out/ugnnum-species_dist.png')
plt.close()

# Correlation of number of genes and associated OGs
fig, ax = plt.subplots()
ax.scatter(h_OG, h_gn)
ax.set_xlabel('Number of Associated OGs')
ax.set_ylabel('Number of Genes')
ax.set_title('Correlation of Numbers of Genes\nand Associated OGs for each Species')

fig.savefig('out/gnnum-OGnum_corr.png')
plt.close()

# Correlation of number of unique genes and associated OGs
fig, ax = plt.subplots()
ax.scatter(h_OG, h_ugn)
ax.set_xlabel('Number of Associated OGs')
ax.set_ylabel('Number of Unique Genes')
ax.set_title('Correlation of Numbers of Unique Genes\nand Associated OGs for each Species')

fig.savefig('out/ugnnum-OGnum_corr.png')
plt.close()

# Distribution of number of species
dist_species = groups['spid'].nunique().value_counts()
spec, spec_count = zip(*dist_species.items())
fig, ax1 = plt.subplots()
ax1.bar(spec, spec_count)
ax1.set_title('Distribution of OGs across Number of Unique Species')
ax1.set_xlabel('Number of Unique Species')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_OGs, mx / num_OGs)
ax2.set_ylabel('Fraction of Total OGs')

fig.tight_layout()
fig.savefig('out/OGnum-spnum_dist.png')
plt.close()

# Distribution of number of genes
dist_seq = groups.size().value_counts()
seq, seq_count = zip(*dist_seq.items())
fig, ax1 = plt.subplots()
ax1.bar(seq, seq_count, width=1, align='center')
ax1.set_title('Distribution of OGs across Number of Genes')
ax1.set_xlabel('Number of Genes')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_OGs, mx / num_OGs)
ax2.set_ylabel('Fraction of Total OGs')

fig.tight_layout()
fig.savefig('out/OGnum-gnnum_dist.png')
plt.close()

# Distribution of number of duplicates
dist_dup = (groups.size() - groups['spid'].nunique()).value_counts()
seq, seq_count = zip(*dist_dup.drop(0, errors='ignore').items())  # Drop 0; ignore error if 0 does not exist
fig, ax1 = plt.subplots()
ax1.bar(seq, seq_count, width=1, align='center')
ax1.set_title('Distribution of OGs across Number of Species Duplicates')
ax1.set_xlabel('Number of Species Duplicates')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / num_OGs, mx / num_OGs)
ax2.set_ylabel('Fraction of Total OGs')

fig.tight_layout()
fig.savefig('out/OGnum-spdup.png')
plt.close()

# Print counts
print('number of OGs:', num_OGs)
print()
print('number of OGs with 10 species:', dist_species[10])
print('fraction of OGs with 10 species:', dist_species[10] / num_OGs)
print()
print('number of OGs with 10 genes:', dist_seq[10])
print('fraction of OGs with 10 genes:', dist_seq[10] / num_OGs)
print()
print('number of OGs with 10 species and 10 genes:', num_u10)
print('fraction of OGs with 10 species and 10 genes:', num_u10 / num_OGs)
print()
print('number of OGs with duplicates:', num_OGs - dist_dup[0])
print('fraction of OGs with duplicates', (num_OGs - dist_dup[0]) / num_OGs)

"""
OUTPUT
number of OGs: 13316

number of OGs with 10 species: 10385
fraction of OGs with 10 species: 0.7798888555121658

number of OGs with 10 genes: 9270
fraction of OGs with 10 genes: 0.6961550015019525

number of OGs with 10 species and 10 genes: 9213
fraction of OGs with 10 species and 10 genes: 0.6918744367677981

number of OGs with duplicates: 1385
fraction of OGs with duplicates 0.10401021327726044

NOTES
These plots are largely based off those in analysis/EggNOGv5_validation/ali_stats/ali_stats.py

DEPENDENCIES
../cluster_xgraph/cluster_ggraph.py
    ../cluster_xgraph/out/gclusters.tsv
../ppid2meta/ppid2meta.py
    ../ppid2meta/out/ppid2meta.tsv
"""