"""Plot various statistics of OGs."""

import os

import matplotlib.pyplot as plt
import pandas as pd

# Load genomes
spids = set()
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spids.add(line.rstrip('\n').split('\t')[0])

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, _ = line.rstrip('\n').split('\t')
        ppid2data[ppid] = gnid, spid

# Load OGs
rows = []
with open('../add_paralogs/out/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, OGid, _, edges = line.rstrip('\n').split('\t')
        ppids = {node for edge in edges.split(',') for node in edge.split(':')}
        for ppid in ppids:
            gnid, spid = ppid2data[ppid]
            rows.append({'component_id': component_id, 'OGid': OGid, 'ppid': ppid, 'gnid': gnid, 'spid': spid})
OGs = pd.DataFrame(rows)

groups = OGs.groupby('OGid')
OGidnum = OGs['OGid'].nunique()
ppidnum = OGs['ppid'].nunique()
gnidnum = OGs['gnid'].nunique()

if not os.path.exists('out/cluster/'):
    os.makedirs('out/cluster/')

# Plots
# 1 DISTRIBUTIONS ACROSS SPECIES
# 1.1 Number of OGs for species
counts_OGid = OGs[['spid', 'OGid']].drop_duplicates()['spid'].value_counts().sort_index()
labels, height = zip(*counts_OGid.items())
xs = list(range(len(labels)))
fig, ax1 = plt.subplots()
ax1.bar(xs, height)
ax1.set_xticks(xs, labels, rotation=60, fontsize=8)
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of associated OGs')

ax2 = ax1.twinx()
ymin, ymax = ax1.get_ylim()
ax2.set_ylim(ymin / OGidnum, ymax / OGidnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/cluster/bar_OGidnum-species.png')
plt.close()

# 1.2 Distribution of proteins across species
counts_ppid = OGs.groupby('spid')['ppid'].nunique().sort_index()
labels, height = zip(*counts_ppid.items())
xs = list(range(len(labels)))
fig, ax1 = plt.subplots()
ax1.bar(xs, height)
ax1.set_xticks(xs, labels, rotation=60, fontsize=8)
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of proteins')

ax2 = ax1.twinx()
ymin, ymax = ax1.get_ylim()
ax2.set_ylim(ymin / ppidnum, ymax / ppidnum)
ax2.set_ylabel('Fraction of total proteins')

fig.tight_layout()
fig.savefig('out/cluster/bar_ppidnum-species.png')
plt.close()

# 1.3 Distribution of genes across species
counts_gnid = OGs.groupby('spid')['gnid'].nunique().sort_index()
labels, height = zip(*counts_gnid.items())
xs = list(range(len(labels)))
fig, ax1 = plt.subplots()
ax1.bar(xs, height)
ax1.set_xticks(xs, labels, rotation=60, fontsize=8)
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of genes')

ax2 = ax1.twinx()
ymin, ymax = ax1.get_ylim()
ax2.set_ylim(ymin / gnidnum, ymax / gnidnum)
ax2.set_ylabel('Fraction of total genes')

fig.tight_layout()
fig.savefig('out/cluster/bar_gnidnum-species.png')
plt.close()

# 1.4 Correlation of number of proteins and associated OGs
fig, ax = plt.subplots()
ax.scatter(counts_OGid, counts_ppid, s=6)
ax.set_xlabel('Number of associated OGs')
ax.set_ylabel('Number of proteins')
ax.set_title('Correlation of numbers of\nproteins and associated OGs for each species')

fig.savefig('out/cluster/scatter_ppidnum-OGidnum.png')
plt.close()

# 1.5 Correlation of number of genes and associated OGs
fig, ax = plt.subplots()
ax.scatter(counts_OGid, counts_gnid, s=6)
ax.set_xlabel('Number of associated OGs')
ax.set_ylabel('Number of genes')
ax.set_title('Correlation of numbers of\ngenes and associated OGs for each species')

fig.savefig('out/cluster/scatter_gnidnum-OGidnum.png')
plt.close()

# 2 NUMBER OF EXCLUSIONS
for i in range(len(spids)-10, len(spids)):
    exclusions = {spid: 0 for spid in spids}
    for _, group in groups:
        included_spids = set(group['spid'])
        if len(included_spids) == i:
            for spid in (spids - included_spids):
                exclusions[spid] += 1
    labels, height = zip(*sorted(exclusions.items()))
    x = list(range(len(labels)))
    fig, ax = plt.subplots()
    ax.bar(x, height)
    ax.set_xticks(xs, labels, rotation=60, fontsize=8)
    ax.set_xlabel('Species')
    ax.set_ylabel('Number of exclusions')
    ax.set_title(f'Number of exclusions in OGs with {i} species')

    fig.tight_layout()
    fig.savefig(f'out/cluster/bar_exclusion-species_{i}.png')
    plt.close()

# 3 DISTRIBUTIONS OF OGS
# 3.1 Distribution of OGs across number of species
counts = groups['spid'].nunique().value_counts()
xs, height = zip(*counts.items())
fig, ax1 = plt.subplots()
ax1.bar(xs, height)
ax1.set_xlabel('Number of species in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
ymin, ymax = ax1.get_ylim()
ax2.set_ylim(ymin / OGidnum, ymax / OGidnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/cluster/hist_OGidnum-spidnum.png')
plt.close()

# 3.2 Distribution of OGs across number of species duplicates
counts = (groups.size() - groups['spid'].nunique()).value_counts()
xs, height = zip(*counts.drop(0, errors='ignore').items())  # Drop 0; ignore error if 0 does not exist
fig, ax1 = plt.subplots()
ax1.bar(xs, height, width=1)
ax1.set_xlabel('Number of species duplicates in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
ymin, ymax = ax1.get_ylim()
ax2.set_ylim(ymin / OGidnum, ymax / OGidnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/cluster/hist_OGidnum-duplicates.png')
plt.close()

# 3.3 Distribution of OGs across number of proteins
counts = groups['ppid'].nunique().value_counts()
xs, height = zip(*counts.items())
fig, ax1 = plt.subplots()
ax1.bar(xs, height, width=1)
ax1.set_xlabel('Number of proteins in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
ymin, ymax = ax1.get_ylim()
ax2.set_ylim(ymin / OGidnum, ymax / OGidnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/cluster/hist_OGidnum-ppidnum.png')
plt.close()

# 3.4 Distribution of OGs across number of genes
counts = groups['gnid'].nunique().value_counts()
xs, height = zip(*counts.items())
fig, ax1 = plt.subplots()
ax1.bar(xs, height, width=1)
ax1.set_xlabel('Number of genes in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
ymin, ymax = ax1.get_ylim()
ax2.set_ylim(ymin / OGidnum, ymax / OGidnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/cluster/hist_OGidnum-gnidnum.png')
plt.close()

# 4 DISTRIBUTIONS ACROSS OGS
# 4.1 Distribution of proteins across number of associated OGs
counts = OGs.groupby('ppid')['OGid'].nunique().value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of OGs associated with protein')
plt.ylabel('Number of proteins')
plt.savefig('out/cluster/hist_ppidnum-OGidnum.png')
plt.yscale('log')
plt.savefig('out/cluster/hist_ppidnum-OGidnum_log.png')
plt.close()

# 4.2 Distribution of genes across number of associated OGs
counts = OGs.groupby('gnid')['OGid'].nunique().value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of OGs associated with gene')
plt.ylabel('Number of genes')
plt.savefig('out/cluster/hist_gnidnum-OGidnum.png')
plt.yscale('log')
plt.savefig('out/cluster/hist_gnidnum-OGidnum_log.png')
plt.close()

# 4.3 Distribution of connected components across number of associated OGs
counts = OGs.groupby('component_id')['OGid'].nunique().value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of OGs in component')
plt.ylabel('Number of components')
plt.savefig('out/cluster/hist_componentnum-OGnum.png')
plt.yscale('log')
plt.savefig('out/cluster/hist_componentnum-OGnum_log.png')
plt.close()

"""
NOTES
DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../config/genomes.tsv
../add_paralogs/add_paralogs.py
    ../add_paralogs/out/clusters.tsv
"""