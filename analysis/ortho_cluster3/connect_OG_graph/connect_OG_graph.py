"""Find connected components of OG graph."""

import os
from src.ortho_cluster.graphs import get_connected_components


def get_sort_tuple(component):
    return len(component), sorted(component)


# Load graph
graph = {}
with open('../OGs2graph/out/OG_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = [adj.split(':')[0] for adj in adjs.split(',')] if adjs else []

# Find connected components
components = get_connected_components(graph)
components = sorted(components, key=get_sort_tuple, reverse=True)

# Write components to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/components.tsv', 'w') as file:
    file.write('GGid\tOGids\n')
    for i, component in enumerate(components):
        GGid = f'{i:04X}'  # Uppercase hex, zero-padded to 4
        nodestring = ','.join(component)
        file.write(f'{GGid}\t{nodestring}\n')

"""
DEPENDENCIES
../OGs2graph/OGs2graph.py
    ../OGs2graph/out/OG_graph.tsv
"""