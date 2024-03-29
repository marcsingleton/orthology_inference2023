"""Find connected components of graph."""

import os
from src.ortho_cluster.graphs import get_connected_components


def get_sort_tuple(component):
    return len(component), sorted(component)


# Load graph
graph = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = [adj.split(':')[0] for adj in adjs.split(',')]

# Find connected components
components = get_connected_components(graph)
components = sorted(filter(lambda x: len(x) > 1, components), key=get_sort_tuple, reverse=True)
# Write clusters to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/components.tsv', 'w') as file:
    file.write('component_id\tppids\n')
    for i, component in enumerate(components):
        component_id = f'{i:04X}'  # Uppercase hex, zero-padded to 4
        nodestring = ','.join(component)
        file.write(f'{component_id}\t{nodestring}\n')

"""
DEPENDENCIES
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""