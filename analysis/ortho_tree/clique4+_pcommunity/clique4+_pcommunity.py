"""Cluster graph by k-clique percolation."""

import os
import signal
from itertools import combinations, product

import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.clique import find_cliques
from networkx.algorithms.components import connected_components
from networkx.algorithms.core import k_core
import src.ortho_cluster.percolate as percolate
from src.ortho_cluster.percolate import edges2nodes


# Output functions
def classify_CC(CCtypes, subOGs):
    subnOGs = []
    for subOG in subOGs:
        subnOG = set()
        for node1, node2 in subOG:
            subnOG.update(sqid2ppids[node1])
            subnOG.update(sqid2ppids[node2])
        subnOGs.append(subnOG)
    if len(subnOGs) == 0:
        CCtypes[0][len(subnOGs)] = CCtypes[0].get(len(subnOGs), 0) + 1  # Component has no OGs
    elif len(subnOGs) == 1:
        if len(subnOGs[0]) == len(CC):
            CCtypes[1][len(subnOGs)] = CCtypes[1].get(len(subnOGs), 0) + 1  # Component and OG are equal
        else:
            CCtypes[2][len(subnOGs)] = CCtypes[2].get(len(subnOGs), 0) + 1  # Component has single OG which is a subset of the component
    elif any([set.intersection(nOG1, nOG2) for nOG1, nOG2 in combinations(subnOGs, 2)]):
        CCtypes[4][len(subnOGs)] = CCtypes[4].get(len(subnOGs), 0) + 1  # Component has multiple non-disjoint OGs
    else:
        CCtypes[3][len(subnOGs)] = CCtypes[3].get(len(subnOGs), 0) + 1  # Component has multiple pairwise disjoint OGs


def save_results(OGs, CCtypes, k):
    # Make plots output directory
    if not os.path.exists(f'out/{k}clique/'):
        os.makedirs(f'out/{k}clique/')  # Recursive folder creation

    # Write OGs to file
    j = 0
    with open(f'out/{k}clique/pclusters.txt', 'w') as file:
        for i, subOGs in enumerate(OGs):
            CCid = hex(i)[2:].zfill(4)
            for OG in sorted(subOGs, key=lambda x: sorted(edges2nodes(x))):  # Ensure consistent ids by keying on sorted node list
                OGid = hex(j)[2:].zfill(4)
                edgestring = '\t'.join([f'{n1},{n2}' for node1, node2 in OG for n1, n2 in product(sqid2ppids[node1], sqid2ppids[node2])])
                file.write(CCid + ':' + OGid + ':' + edgestring + '\n')
                j += 1

    # Plots
    plt.bar(CCtypes[0].keys(), CCtypes[0].values(), label='Type 0')
    plt.bar(CCtypes[1].keys(), CCtypes[1].values(), label='Type 1')
    plt.bar(CCtypes[2].keys(), CCtypes[2].values(), bottom=CCtypes[1].get(1, 0), label='Type 2')
    plt.bar(CCtypes[3].keys(), CCtypes[3].values(), label='Type 3')
    plt.bar(CCtypes[4].keys(), CCtypes[4].values(), bottom=[CCtypes[3].get(key, 0) for key in CCtypes[4]], label='Type 4')
    plt.xlabel('Number of OGs in connected component')
    plt.ylabel('Number of connected components')
    plt.title('Distribution of connected components across number of OGs')
    plt.legend()
    plt.savefig(f'out/{k}clique/connectnum-OGnum_type_dist1-1.png')
    plt.xlim((-1, 17))  # Adjust axis to truncate outliers
    plt.savefig(f'out/{k}clique/connectnum-OGnum_type_dist1-2.png')
    plt.close()

    plt.bar(CCtypes[3].keys(), CCtypes[3].values(), label='Type 3', color='C3')
    plt.bar(CCtypes[4].keys(), CCtypes[4].values(), bottom=[CCtypes[3].get(key, 0) for key in CCtypes[4]], label='Type 4', color='C4')
    plt.xlabel('Number of OGs in connected component')
    plt.ylabel('Number of connected components')
    plt.title('Distribution of connected components across number of OGs')
    plt.legend()
    plt.savefig(f'out/{k}clique//connectnum-OGnum_type_dist2-1.png')
    plt.xlim((-1, 17))  # Adjust axis to truncate outliers
    plt.savefig(f'out/{k}clique/connectnum-OGnum_type_dist2-2.png')
    plt.close()

    plt.pie([sum(CCtype.values()) for CCtype in CCtypes], labels=[f'Type {i}' for i in range(len(CCtypes))])
    plt.title('Connected components by type')
    plt.savefig(f'out/{k}clique/type_pie.png')
    plt.close()

    print()
    print(f'{k}-CLIQUE')
    for i, CCtype in enumerate(CCtypes):
        print(f'Type {i}:', sum(CCtype.values()))


# Load seq metadata
sqid2ppids = {}
constituents = set()  # Constituent ppids
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, _, _, sqid = line.split()
        try:
            sqid2ppids[sqid].append(ppid)
        except KeyError:
            sqid2ppids[sqid] = [ppid]
        if ppid != sqid:
            constituents.add(ppid)

# Load pgraph
graph = {}
with open('../hits2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = set([adj.split(':')[0] for adj in adjs.split(',')])

# Load connected components
CCs = []
with open('../connect_pgraph/out/pconnect.txt') as file:
    for line in file:
        CCid, nodes = line.rstrip().split(':')
        CCs.append((CCid, set(nodes.split(','))))

ks = list(range(4, 7))
OGs_ks = {k: [] for k in ks}
CCtypes_ks = {k: [{} for _ in range(5)] for k in ks}
for CCid, CC in CCs:
    # Create graph
    G = nx.Graph()
    for node in CC:
        if node in constituents:
            continue
        G.add_node(node)
        for adj in filter(lambda x: x not in constituents, graph[node]):
            G.add_edge(node, adj)

    # Handle cliques
    try:
        signal.signal(signal.SIGALRM, percolate.clique_handler)
        if CCid in ['03ae', '03b0', '03b2', '03b5', '03b7', '0893']:  # Skip histones and other complex OGs
            raise percolate.CliqueError
        signal.alarm(90)
        cliques = list(find_cliques(G))
        signal.alarm(0)
    except percolate.CliqueError:
        print(f'CliqueError: {CCid}')
        for k in ks:
            subOGs = set()
            core = k_core(G, k)
            for component in connected_components(core):
                subOGs.add(frozenset([frozenset(edge) for edge in core.edges(component)]))
            OGs_ks[k].append(subOGs)
            classify_CC(CCtypes_ks[k], subOGs)
        continue  # Continue to next OG

    # Handle percolation
    for k in ks:
        try:
            signal.signal(signal.SIGALRM, percolate.percolate_handler)
            signal.alarm(90)
            subOGs = list(percolate.k_clique_communities_progressive(G, k, cliques))
            signal.alarm(0)
        except percolate.PercolateError:
            print(f'PercolateError: ({k}, {CCid})')
            subOGs = set()
            core = k_core(G, k)
            for component in connected_components(core):
                subOGs.add(frozenset([frozenset(edge) for edge in core.edges(component)]))
        OGs_ks[k].append(subOGs)
        classify_CC(CCtypes_ks[k], subOGs)

for k in ks:
    save_results(OGs_ks[k], CCtypes_ks[k], k)

"""
OUTPUT
CliqueError: 03ae
CliqueError: 03b0
CliqueError: 03b2
CliqueError: 03b5
CliqueError: 03b7
CliqueError: 0893

4-CLIQUE
Type 0: 3214
Type 1: 10982
Type 2: 3858
Type 3: 646
Type 4: 1908

5-CLIQUE
Type 0: 4443
Type 1: 9647
Type 2: 4069
Type 3: 663
Type 4: 1786

6-CLIQUE
Type 0: 5279
Type 1: 8759
Type 2: 4223
Type 3: 720
Type 4: 1627

DEPENDENCIES
../connect_pgraph/connect_pgraph.py
    ../connect_pgraph/out/pconnect.txt
../hits2pgraph/hits2pgraph.py
    ../hits2pgraph/out/pgraph.tsv
"""