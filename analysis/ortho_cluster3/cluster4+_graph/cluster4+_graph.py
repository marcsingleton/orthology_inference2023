"""Cluster graph by k-clique percolation."""

import os
import signal
from itertools import combinations

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import src.ortho_cluster.percolate as percolate
from networkx.algorithms.clique import find_cliques
from networkx.algorithms.components import connected_components
from networkx.algorithms.core import k_core


def get_component_type(OGs, component):
    node_sets = [{node for edge in OG for node in edge} for OG in OGs]
    if len(node_sets) == 0:
        return 0  # Component has no OGs
    elif len(node_sets) == 1:
        if len(node_sets[0]) == len(component):
            return 1  # Component and OG are equal
        else:
            return 2  # Component has single OG which is a subset of the component
    elif any([set.intersection(set1, set2) for set1, set2 in combinations(node_sets, 2)]):
        return 4  # Component has multiple non-disjoint OGs
    else:
        return 3  # Component has multiple pairwise disjoint OGs


def get_sort_tuple(OG_record):
    OG = OG_record[1]
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
OG_records_ks = [[] for _ in ks]
rows_ks = [[] for _ in ks]
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
        for k, OG_records_k, rows_k in zip(ks, OG_records_ks, rows_ks):
            core = k_core(subgraph, k)
            OGs = [core.edges(core_component) for core_component in connected_components(core)]
            OG_records_k.extend([(component_id, OG, f'{k}core') for OG in OGs])
            component_type = get_component_type(OGs, component)
            rows_k.append({'component_id': component_id, 'component_type': component_type, 'OGnum': len(OGs)})
        continue  # Continue to next OG

    # Handle percolation
    for k, OG_records_k, rows_k in zip(ks, OG_records_ks, rows_ks):
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
        OG_records_k.extend([(component_id, OG, algorithm) for OG in OGs])
        component_type = get_component_type(OGs, component)
        rows_k.append({'component_id': component_id, 'component_type': component_type, 'OGnum': len(OGs)})

for k, OG_records_k, rows_k in zip(ks, OG_records_ks, rows_ks):
    if not os.path.exists(f'out/{k}clique/'):
        os.makedirs(f'out/{k}clique/')

    # Write OGs to file
    j = 0
    with open(f'out/{k}clique/clusters.tsv', 'w') as file:
        file.write('component_id\tOGid\talgorithm\tedges\n')
        for component_id, OG, algorithm in sorted(OG_records_k, key=get_sort_tuple, reverse=True):
            OGid = hex(j)[2:].zfill(4).upper()
            edgestring = ','.join([f'{node1}:{node2}' for node1, node2 in OG])
            file.write(f'{component_id}\t{OGid}\t{algorithm}\t{edgestring}\n')
            j += 1

    # Plots
    df = pd.DataFrame(rows_k)
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
    plt.legend()
    plt.savefig(f'out/{k}clique/bar_connectnum-OGnum_type_dist1-1.png')
    plt.xlim((-1, 17))  # Adjust axis to truncate outliers
    plt.savefig(f'out/{k}clique/bar_connectnum-OGnum_type_dist1-2.png')
    plt.close()

    plt.bar(component_types[3].index, component_types[3].values, label='Type 3', color='C3')
    plt.bar(component_types[4].index, component_types[4].values,
            bottom=[component_types[3].get(index, 0) for index in component_types[4].index], label='Type 4', color='C4')
    plt.xlabel('Number of OGs in connected component')
    plt.ylabel('Number of connected components')
    plt.legend()
    plt.savefig(f'out/{k}clique/bar_connectnum-OGnum_type_dist2-1.png')
    plt.xlim((-1, 17))  # Adjust axis to truncate outliers
    plt.savefig(f'out/{k}clique/bar_connectnum-OGnum_type_dist2-2.png')
    plt.close()

    counts = [component_type.sum() for component_type in component_types]
    labels = [f'Type {i}\n{count}' for i, count in zip(range(len(component_types)), counts)]
    plt.pie(counts, labels=labels, labeldistance=1.3, textprops={'ha': 'center'})
    plt.title('Connected components by type')
    plt.savefig('out/pie_component_type.png')
    plt.close()

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