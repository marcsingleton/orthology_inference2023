"""Extract clusters from connected components of graph."""

import os
from src.ortho_cluster.graphs import get_triangle_clusters


def get_sort_tuple(record):
    OG = record[1]
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

# Cluster by triangle criterion
records = []
for component_id, component in components:
    subgraph = {node: graph[node] for node in component}

    OGs = get_triangle_clusters(subgraph)
    records.extend([(component_id, OG, '3clique') for OG in OGs])

# Write OGs to file
if not os.path.exists('out/'):
    os.mkdir('out/')

j = 0
with open('out/clusters.tsv', 'w') as file:
    file.write('component_id\tOGid\talgorithm\tedges\n')
    for component_id, OG, algorithm in sorted(records, key=get_sort_tuple, reverse=True):
        OGid = hex(j)[2:].zfill(4).upper()
        edgestring = ','.join([f'{node1}:{node2}' for node1, node2 in OG])
        file.write(f'{component_id}\t{OGid}\t{algorithm}\t{edgestring}\n')
        j += 1

"""
DEPENDENCIES
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""