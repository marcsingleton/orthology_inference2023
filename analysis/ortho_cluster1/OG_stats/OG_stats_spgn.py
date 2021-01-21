"""Plot various statistics of OGs relating to their counts of genes and species."""

import matplotlib.pyplot as plt
import os
import pandas as pd

# Load seq metadata
gnid2spid = {}
with open('../ppid2meta/out/seq_meta.tsv') as file:
    for line in file:
        _, gnid, spid = line.split()
        gnid2spid[gnid] = spid

# Load OGs
rows = []
with open('../subcluster_xgraph/out/ggraph/gclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for gnid in gnids:
            rows.append({'CCid': CCid, 'OGid': OGid, 'gnid': gnid, 'spid': gnid2spid[gnid]})
OGs = pd.DataFrame(rows)

groups = OGs.groupby('OGid')
OGnum = OGs['OGid'].nunique()
ugnnum = OGs['gnid'].nunique()
u10num = len(groups.filter(lambda x: len(x) == 10 and x['spid'].nunique() == 10)) // 10

# Make output directory
if not os.path.exists('out/spgn/'):
    os.makedirs('out/spgn/')  # Recursive folder creation

# Plots
# Number of associated OGs for species
spid_OGs = {}
for _, group in groups:
    for spid in group['spid'].unique():
        spid_OGs[spid] = spid_OGs.get(spid, 0) + 1

labels, h_OG = zip(*sorted(spid_OGs.items(), key=lambda i: i[0]))
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_OG, tick_label=labels)
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of associated OGs')
ax1.set_title('Number of associated OGs for each species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/spgn/bar_OGnum-species.png')
plt.close()

# Distribution of genes across species
spid_gns = OGs['spid'].value_counts().sort_index()
labels, h_gn = zip(*spid_gns.items())
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_gn, tick_label=labels)
ax1.set_xlabel('Associated species')
ax1.set_ylabel('Number of genes')
ax1.set_title('Distribution of genes across associated species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / len(OGs), mx / len(OGs))
ax2.set_ylabel('Fraction of total genes')

fig.tight_layout()
fig.savefig('out/spgn/bar_gnnum-species.png')
plt.close()

# Distribution of unique genes across species
spid_ugns = OGs.groupby('spid')['gnid'].nunique().sort_index()
labels, h_ugn = zip(*spid_ugns.items())
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_ugn, tick_label=labels)
ax1.set_xlabel('Associated species')
ax1.set_ylabel('Number of unique genes')
ax1.set_title('Distribution of unique genes across associated species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / ugnnum, mx / ugnnum)
ax2.set_ylabel('Fraction of total unique genes')

fig.tight_layout()
fig.savefig('out/spgn/bar_ugnnum-species.png')
plt.close()

# Number of exclusions for each species
spids = set(OGs['spid'].drop_duplicates())
spid_counts = {spid: 0 for spid in sorted(spids)}
for spid, in [spids - set(group['spid'].drop_duplicates()) for _, group in groups if group['spid'].nunique() == 9]:
    spid_counts[spid] += 1
labels, h = zip(*spid_counts.items())
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h, tick_label=labels)
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of OGs')
ax1.set_title('Number of exclusions for each species in OGs with 9 species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / sum(h), mx / sum(h))
ax2.set_ylabel('Fraction of total OGs with 9 species')

fig.tight_layout()
fig.savefig('out/spgn/bar_OGexclusion-species.png')
plt.close()

# Correlation of number of genes and associated OGs
fig, ax = plt.subplots()
ax.scatter(h_OG, h_gn)
ax.set_xlabel('Number of associated OGs')
ax.set_ylabel('Number of associated genes')
ax.set_title('Correlation of numbers of associated\ngenes and OGs for each species')

fig.savefig('out/spgn/scatter_gnnum-OGnum.png')
plt.close()

# Correlation of number of unique genes and associated OGs
fig, ax = plt.subplots()
ax.scatter(h_OG, h_ugn)
ax.set_xlabel('Number of associated OGs')
ax.set_ylabel('Number of associated unique genes')
ax.set_title('Correlation of numbers of associated\n unique genes and OGs for each species')

fig.savefig('out/spgn/scatter_ugnnum-OGnum.png')
plt.close()

# Distribution of OGs across number of unique species
dist_species = groups['spid'].nunique().value_counts()
spec, spec_count = zip(*dist_species.items())
fig, ax1 = plt.subplots()
ax1.bar(spec, spec_count)
ax1.set_title('Distribution of OGs across number of unique species in OG')
ax1.set_xlabel('Number of unique species in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/spgn/hist_OGnum-spnum.png')
plt.close()

# Distribution of OGs across number of genes
dist_seq = groups.size().value_counts()
seq, seq_count = zip(*dist_seq.items())
fig, ax1 = plt.subplots()
ax1.bar(seq, seq_count, width=1)
ax1.set_title('Distribution of OGs across number of genes in OG')
ax1.set_xlabel('Number of genes in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/spgn/hist_OGnum-gnnum.png')
plt.close()

# Distribution of OGs across number of species duplicates
dist_dup = (groups.size() - groups['spid'].nunique()).value_counts()
seq, seq_count = zip(*dist_dup.drop(0, errors='ignore').items())  # Drop 0; ignore error if 0 does not exist
fig, ax1 = plt.subplots()
ax1.bar(seq, seq_count, width=1)
ax1.set_title('Distribution of OGs across number of species duplicates in OG')
ax1.set_xlabel('Number of species duplicates in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/spgn/hist_OGnum-spdup.png')
plt.close()

# Print counts
print('number of OGs:', OGnum)
print()
print('number of OGs with 10 species:', dist_species[10])
print('fraction of OGs with 10 species:', dist_species[10] / OGnum)
print()
print('number of OGs with 10 genes:', dist_seq[10])
print('fraction of OGs with 10 genes:', dist_seq[10] / OGnum)
print()
print('number of OGs with 10 species and 10 genes:', u10num)
print('fraction of OGs with 10 species and 10 genes:', u10num / OGnum)
print()
print('number of OGs with duplicates:', OGnum - dist_dup[0])
print('fraction of OGs with duplicates', (OGnum - dist_dup[0]) / OGnum)

"""
OUTPUT
number of OGs: 13305

number of OGs with 10 species: 10385
fraction of OGs with 10 species: 0.7805336339721909

number of OGs with 10 genes: 8886
fraction of OGs with 10 genes: 0.667869222096956

number of OGs with 10 species and 10 genes: 8797
fraction of OGs with 10 species and 10 genes: 0.6611800075159714

number of OGs with duplicates: 1924
fraction of OGs with duplicates 0.14460729049229612

NOTES
These plots are largely based off those in analysis/EggNOGv5_validation/ali_stats/ali_stats.py

DEPENDENCIES
../seq_meta/seq_meta.py
    ../seq_meta/out/seq_meta.tsv
../subcluster_xgraph/subcluster_ggraph.py
    ../subcluster_xgraph/out/gclusters.tsv
"""