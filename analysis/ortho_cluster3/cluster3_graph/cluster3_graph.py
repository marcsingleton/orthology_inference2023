"""Extract clusters from connected components of graph."""

import os
from itertools import combinations

import matplotlib.pyplot as plt
import pandas as pd
from src.ortho_cluster.graphs import get_triangle_clusters


def get_sort_tuple(OG_record):
    OG = OG_record[1]
    nodes = sorted([node for edge in OG for node in edge])
    return len(nodes), nodes


# Load graph
graph = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = {adj.split(':')[0] for adj in adjs.split(',')}

# Load connected components
components = []
with open('../connect_hit_graph/out/components.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        component_id, nodes = line.rstrip('\n').split('\t')
        components.append((component_id, set(nodes.split(','))))

OG_records, rows = [], []
for component_id, component in components:
    subgraph = {node: graph[node] for node in component}

    # Cluster by triangle criterion
    OGs = get_triangle_clusters(subgraph)
    OG_records.extend([(component_id, OG, '3clique') for OG in OGs])

    # Classify component
    node_sets = [{node for edge in OG for node in edge} for OG in OGs]
    if len(node_sets) == 0:
        component_type = 0  # Component has no OGs
    elif len(node_sets) == 1:
        if len(node_sets[0]) == len(component):
            component_type = 1  # Component and OG are equal
        else:
            component_type = 2  # Component has single OG which is a subset of the component
    elif any([set.intersection(set1, set2) for set1, set2 in combinations(node_sets, 2)]):
        component_type = 4  # Component has multiple non-disjoint OGs
    else:
        component_type = 3  # Component has multiple pairwise disjoint OGs

    rows.append({'component_id': component_id, 'component_type': component_type, 'OGnum': len(OGs)})

if not os.path.exists('out/'):
    os.mkdir('out/')

# Write OGs to file
j = 0
with open('out/clusters.tsv', 'w') as file:
    file.write('component_id\tOGid\talgorithm\tedges\n')
    for component_id, OG, algorithm in sorted(OG_records, key=get_sort_tuple, reverse=True):
        OGid = hex(j)[2:].zfill(4).upper()
        edgestring = ','.join([f'{node1}:{node2}' for node1, node2 in OG])
        file.write(f'{component_id}\t{OGid}\t{algorithm}\t{edgestring}\n')
        j += 1

# Plots
df = pd.DataFrame(rows)
component_types = [df.loc[df['component_type'] == i, 'OGnum'].value_counts() for i in range(5)]

plt.bar(component_types[0].index, component_types[0].values, label='Type 0')
plt.bar(component_types[1].index, component_types[1].values, label='Type 1')
plt.bar(component_types[2].index, component_types[2].values,
        bottom=component_types[1].get(1, 0), label='Type 2')
plt.bar(component_types[3].index, component_types[3].values, label='Type 3')
plt.bar(component_types[4].index, component_types[4].values,
        bottom=[component_types[3].get(index, 0) for index in component_types[4].index], label='Type 4')
plt.xlabel('Number of OGs in connected component')
plt.ylabel('Number of connected components')
plt.title('Distribution of connected components across number of OGs')
plt.legend()
plt.savefig('out/bar_connectnum-OGnum_type_dist1-1.png')
plt.xlim((-1, 17))  # Adjust axis to truncate outliers
plt.savefig('out/bar_connectnum-OGnum_type_dist1-2.png')
plt.close()

plt.bar(component_types[3].index, component_types[3].values, label='Type 3', color='C3')
plt.bar(component_types[4].index, component_types[4].values,
        bottom=[component_types[3].get(index, 0) for index in component_types[4].index], label='Type 4', color='C4')
plt.xlabel('Number of OGs in connected component')
plt.ylabel('Number of connected components')
plt.title('Distribution of connected components across number of OGs')
plt.legend()
plt.savefig('out/bar_connectnum-OGnum_type_dist2-1.png')
plt.xlim((-1, 17))  # Adjust axis to truncate outliers
plt.savefig('out/bar_connectnum-OGnum_type_dist2-2.png')
plt.close()

plt.pie([sum(component_type.values) for component_type in component_types],
        labels=[f'Type {i}' for i in range(len(component_types))])
plt.title('Connected components by type')
plt.savefig('out/pie_component_type.png')
plt.close()

for i, component_type in enumerate(component_types):
    print(f'Type {i}:', sum(component_type.values))

"""
OUTPUT
Type 0: 1199
Type 1: 14031
Type 2: 3117
Type 3: 514
Type 4: 1918

DEPENDENCIES
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""