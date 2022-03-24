"""Plot various statistics of OGs relating to counts of connected components and OGs."""

import os

import matplotlib.pyplot as plt
import pandas as pd

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, _ = line.split()
        ppid2meta[ppid] = gnid, spid

# Load connected components
rows = []
with open('../connect_hit_graph/out/components.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, nodes = line.rstrip().split('\t')
        for ppid in nodes.split(','):
            rows.append({'component_id': component_id, 'ppid': ppid})
components = pd.DataFrame(rows)

# Load OGs
rows = []
with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, OGid, _, edges = line.rstrip().split('\t')
        ppids = {node for edge in edges.split(',') for node in edge.split(':')}
        for ppid in ppids:
            gnid, spid = ppid2meta[ppid]
            rows.append({'component_id': component_id, 'OGid': OGid, 'ppid': ppid, 'gnid': gnid, 'spid': spid})
OGs = pd.DataFrame(rows)

# Make output directory
if not os.path.exists('out/CCOG/'):
    os.makedirs('out/CCOG/')  # Recursive folder creation

# Plots
# Distribution of proteins across number of associated OGs
groups = OGs.groupby('ppid')
counts = groups['OGid'].nunique().value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of OGs associated with protein')
plt.ylabel('Number of proteins')
plt.savefig('out/CCOG/hist_ppidnum-OGidnum.png')
plt.yscale('log')
plt.savefig('out/CCOG/hist_ppidnum-OGidnum_log.png')
plt.close()

# Distribution of genes across number of associated OGs
groups = OGs.groupby('gnid')
counts = groups['OGid'].nunique().value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of OGs associated with gene')
plt.ylabel('Number of genes')
plt.savefig('out/CCOG/hist_gnidnum-OGidnum.png')
plt.yscale('log')
plt.savefig('out/CCOG/hist_gnidnum-OGidnum_log.png')
plt.close()

# Distribution of connected components across number of associated OGs
groups = OGs.groupby('component_id')
counts = groups['OGid'].nunique().value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of OGs in component')
plt.ylabel('Number of components')
plt.savefig('out/CCOG/hist_componentnum-OGnum.png')
plt.yscale('log')
plt.savefig('out/CCOG/hist_componentnum-OGnum_log.png')
plt.close()

"""
DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../cluster4+_graph/cluster4+_graph.py
    ../cluster4+_graph/out/4clique/clusters.tsv
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
"""