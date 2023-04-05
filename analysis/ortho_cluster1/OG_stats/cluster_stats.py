"""Plot various statistics of OGs."""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

# Load genomes
spids = set()
with open('../config/genomes.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        spids.add(fields['spid'])

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppid2data[fields['ppid']] = (fields['gnid'], fields['spid'])

# Load OGs
rows = []
with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        component_id, OGid = fields['component_id'], fields['OGid']
        ppids = {node for edge in fields['edges'].split(',') for node in edge.split(':')}
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
ax.scatter(counts_OGid, counts_ppid, s=10)
ax.set_xlabel('Number of associated OGs')
ax.set_ylabel('Number of proteins')
ax.set_title('Correlation of numbers of\nproteins and associated OGs for each species')

fig.savefig('out/cluster/scatter_ppidnum-OGidnum.png')
plt.close()

# 1.5 Correlation of number of genes and associated OGs
fig, ax = plt.subplots()
ax.scatter(counts_OGid, counts_gnid, s=10)
ax.set_xlabel('Number of associated OGs')
ax.set_ylabel('Number of genes')
ax.set_title('Correlation of numbers of\ngenes and associated OGs for each species')

fig.savefig('out/cluster/scatter_gnidnum-OGidnum.png')
plt.close()

# 2 OG MEMBERSHIP PLOTS
# 2.1 Number of exclusions
exclusions = {spid: 0 for spid in spids}
for _, group in groups:
    included_spids = set(group['spid'])
    if len(included_spids) == len(spids) - 1:
        for spid in (spids - included_spids):
            exclusions[spid] += 1
labels, height = zip(*sorted(exclusions.items()))
xs = list(range(len(labels)))
fig, ax = plt.subplots()
ax.bar(xs, height)
ax.set_xticks(xs, labels, rotation=60, fontsize=8)
ax.set_xlabel('Species')
ax.set_ylabel('Number of exclusions')

fig.tight_layout()
fig.savefig(f'out/cluster/bar_exclusion-species.png')
plt.close()

# 2.2 PCA of OG membership
spid2idx = {spid: i for i, spid in enumerate(spids)}
X = np.zeros((len(spids), len(groups)))
for i, (_, group) in enumerate(groups):
    for spid in group['spid']:
        X[spid2idx[spid], i] = 1
pca = PCA().fit(X)
Y = pca.transform(X)

clusters = [('dalb', 'dari', 'dinn', 'dgri', 'dnov', 'dvir', 'dhyd', 'dmoj', 'dnav'),
            ('dwil',),
            ('dper', 'dpse', 'dobs', 'dgua', 'dsob'),
            ('dana', 'dbip', 'dkik', 'dser'),
            ('dfic', 'dele', 'drho', 'dtak', 'dbia', 'dspu', 'dsuz', 'deug'),
            ('dere', 'dtei', 'dsan', 'dyak', 'dmel', 'dmau', 'dsec', 'dsim')]
fig, ax = plt.subplots()
for cluster in clusters:
    idxs = [spid2idx[spid] for spid in cluster]
    ax.scatter(Y[idxs, 0], Y[idxs, 1], label='\n'.join(cluster), alpha=0.9, edgecolors='none')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
fig.subplots_adjust(right=0.85)
fig.savefig('out/cluster/scatter_membership.png')
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
DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../config/genomes.tsv
../cluster4+_graph/cluster4+_graph.py
    ../cluster4+_graph/out/4clique/clusters.tsv
"""