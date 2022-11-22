"""Plot various statistics of components."""

from math import ceil

import matplotlib.pyplot as plt
import pandas as pd

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppid2data[fields['ppid']] = (fields['gnid'], fields['spid'])

graph = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = [adj.split(':') for adj in adjs.split(',')]

# Load connected components
components = {}
with open('out/components.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        components[fields['component_id']] = set(fields['ppids'].split(','))

# Calculate component stats
rows = []
for component_id, component in components.items():
    ppidnum = len(component)
    gnidnum = len(set([ppid2data[ppid][0] for ppid in component]))
    spidnum = len(set([ppid2data[ppid][1] for ppid in component]))
    rows.append({'component_id': component_id, 'ppidnum': ppidnum, 'gnidnum': gnidnum, 'spidnum': spidnum})
df = pd.DataFrame(rows)
df.to_csv('out/stats.tsv', sep='\t', index=False)

degrees = sorted([len(adjs) for adjs in graph.values()])
idx = ceil(len(degrees) * 0.99)
counts = {}
for degree in degrees[:idx]:
    counts[degree] = counts.get(degree, 0) + 1

plt.hexbin(df['ppidnum'], df['gnidnum'], bins='log', gridsize=50, mincnt=1, linewidth=0)
plt.xlabel('Number of proteins in component')
plt.ylabel('Number of unique genes in component')
plt.colorbar()
plt.savefig('out/hexbin_gnidnum-ppidnum.png')
plt.close()

plt.hexbin(df['ppidnum'], df['gnidnum'], bins='log', gridsize=50, mincnt=1, linewidth=0)
plt.xlabel('Number of proteins in component')
plt.ylabel('Number of unique genes in component')
plt.colorbar()
plt.savefig('out/hexbin_spidnum-ppidnum.png')
plt.close()

plt.bar(counts.keys(), counts.values(), width=1)
plt.xlabel('Degree of node')
plt.ylabel('Number of nodes')
plt.savefig('out/hist_nodenum-degree.png')
plt.close()

"""
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
./connect_hit_graph.py
    ./out/components.tsv
"""