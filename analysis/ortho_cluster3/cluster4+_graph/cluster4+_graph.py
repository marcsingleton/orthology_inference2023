"""Cluster graph by k-clique percolation."""

import os
import signal

import networkx as nx
import src.ortho_cluster.percolate as percolate
from networkx.algorithms.clique import find_cliques
from networkx.algorithms.components import connected_components
from networkx.algorithms.core import k_core


def get_sort_tuple(record):
    OG = record[1]
    nodes = sorted([node for edge in OG for node in edge])
    return len(nodes), nodes


clique_timeout = 90
percolate_timeout = 90
k_core_components = ['0008', '0009', '000E', '001D']  # Histones and other complex OGs which directly use k_core

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

ks = list(range(4, 7))
records_ks = [[] for _ in ks]
for component_id, component in components:
    # Create subgraph
    subgraph = nx.Graph()
    for node in component:
        subgraph.add_node(node)
        for adj in graph[node]:
            subgraph.add_edge(node, adj)

    # Handle cliques
    try:
        if component_id in k_core_components:
            raise percolate.CliqueError
        signal.signal(signal.SIGALRM, percolate.clique_handler)
        signal.alarm(clique_timeout)
        cliques = list(find_cliques(subgraph))
        signal.alarm(0)
    except percolate.CliqueError:
        print(f'CliqueError: {component_id}')
        for k, records_k in zip(ks, records_ks):
            core = k_core(subgraph, k)
            OGs = [core.edges(core_component) for core_component in connected_components(core)]
            records_k.extend([(component_id, OG, f'{k}core') for OG in OGs])
        continue  # Continue to next OG

    # Handle percolation
    for k, records_k in zip(ks, records_ks):
        try:
            signal.signal(signal.SIGALRM, percolate.percolate_handler)
            signal.alarm(percolate_timeout)
            OGs = list(percolate.k_clique_communities_progressive(subgraph, k, cliques))
            signal.alarm(0)
            algorithm = f'{k}clique'
        except percolate.PercolateError:
            print(f'PercolateError: ({k}, {component_id})')
            core = k_core(subgraph, k)
            OGs = [core.edges(core_component) for core_component in connected_components(core)]
            algorithm = f'{k}core'
        records_k.extend([(component_id, OG, algorithm) for OG in OGs])

# Write OGs to file
for k, records_k in zip(ks, records_ks):
    if not os.path.exists(f'out/{k}clique/'):
        os.makedirs(f'out/{k}clique/')

    j = 0
    with open(f'out/{k}clique/clusters.tsv', 'w') as file:
        file.write('component_id\tOGid\talgorithm\tedges\n')
        for component_id, OG, algorithm in sorted(records_k, key=get_sort_tuple, reverse=True):
            OGid = hex(j)[2:].zfill(4).upper()
            edgestring = ','.join([f'{node1}:{node2}' for node1, node2 in OG])
            file.write(f'{component_id}\t{OGid}\t{algorithm}\t{edgestring}\n')
            j += 1

"""
OUTPUT
CliqueError: 0008
CliqueError: 0009
CliqueError: 000E
CliqueError: 001D
CliqueError: 00A2
PercolateError: (4, 02D0)
PercolateError: (5, 02D0)
PercolateError: (6, 02D0)
CliqueError: 0486

DEPENDENCIES
../connect_hit_graph/connect_hit_graph.py
    ../connect_hit_graph/out/components.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""