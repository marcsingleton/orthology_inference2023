"""Find connected components of graph."""

import os
from src.ortho_cluster.graphs import get_connected_components

# Load graph
graph = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = [adj.split(':')[0] for adj in adjs.split(',')]

# Find connected components
components = get_connected_components(graph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/components.tsv', 'w') as file:
    file.write('component_id\tppids\n')
    for i, component in enumerate(filter(lambda x: len(x) > 1, components)):
        component_id = hex(i)[2:].zfill(4).upper()
        file.write(component_id + '\t' + ','.join(component) + '\n')

"""
DEPENDENCIES
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""