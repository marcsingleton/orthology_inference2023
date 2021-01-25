"""Cluster polypeptide graphs within OGs by k-clique percolation."""

import os
import signal
from itertools import combinations

import matplotlib.pyplot as plt
import networkx as nx
import src.ortho_cluster.percolate as percolate
from networkx.algorithms.clique import find_cliques
from networkx.algorithms.components import connected_components
from networkx.algorithms.core import k_core
from src.ortho_cluster.percolate import edges2nodes


# Output functions
def classify_OG(OGtypes, subpOGs):
    subnOGs = [edges2nodes(subpOG) for subpOG in subpOGs]
    if len(subnOGs) == 0:
        OGtypes[0][len(subnOGs)] = OGtypes[0].get(len(subnOGs), 0) + 1  # OG has no pOGs
    elif len(subnOGs) == 1:
        if len(set([ppid2gnid[ppid] for ppid in subnOGs[0]])) == len(OG):
            OGtypes[1][len(subnOGs)] = OGtypes[1].get(len(subnOGs), 0) + 1  # OG and pOG are equal
        else:
            OGtypes[2][len(subnOGs)] = OGtypes[2].get(len(subnOGs), 0) + 1  # OG has single pOG which is a subset of the OG
    elif any([set.intersection(nOG1, nOG2) for nOG1, nOG2 in combinations(subnOGs, 2)]):
        OGtypes[4][len(subnOGs)] = OGtypes[4].get(len(subnOGs), 0) + 1  # OG has multiple non-disjoint pOGs
    else:
        OGtypes[3][len(subnOGs)] = OGtypes[3].get(len(subnOGs), 0) + 1  # OG has multiple pairwise disjoint pOGs


def save_results(pOGs, OGtypes, k):
    # Make plots output directory
    if not os.path.exists(f'out/{k}clique/'):
        os.makedirs(f'out/{k}clique/')  # Recursive folder creation

    # Write pOGs to file
    j = 0
    with open(f'out/{k}clique/pclusters.txt', 'w') as file:
        for i, subpOGs in enumerate(pOGs):
            OGid = hex(i)[2:].zfill(4)
            for pOG in sorted(subpOGs, key=lambda x: sorted(edges2nodes(x))):  # Ensure consistent ids by keying on sorted node list
                pOGid = hex(j)[2:].zfill(4)
                file.write(OGid + ':' + pOGid + ':' + '\t'.join([f'{node1},{node2}' for node1, node2 in pOG]) + '\n')
                j += 1

    # Plots
    plt.bar(OGtypes[0].keys(), OGtypes[0].values(), label='Type 0')
    plt.bar(OGtypes[1].keys(), OGtypes[1].values(), label='Type 1')
    plt.bar(OGtypes[2].keys(), OGtypes[2].values(), bottom=OGtypes[1][1], label='Type 2')
    plt.bar(OGtypes[3].keys(), OGtypes[3].values(), label='Type 3')
    plt.bar(OGtypes[4].keys(), OGtypes[4].values(), bottom=[OGtypes[3].get(key, 0) for key in OGtypes[4]], label='Type 4')
    plt.xlabel('Number of pOGs in OG')
    plt.ylabel('Number of OGs')
    plt.title('Distribution of OGs across number of pOGs')
    plt.legend()
    plt.savefig(f'out/{k}clique/OGnum-pOGnum_type_dist1-1.png')
    plt.xlim((-1, 17))  # Adjust axis to truncate outliers
    plt.savefig(f'out/{k}clique/OGnum-pOGnum_type_dist1-2.png')
    plt.close()

    plt.bar(OGtypes[3].keys(), OGtypes[3].values(), label='Type 3', color='C3')
    plt.bar(OGtypes[4].keys(), OGtypes[4].values(), bottom=[OGtypes[3].get(key, 0) for key in OGtypes[4]], label='Type 4', color='C4')
    plt.xlabel('Number of pOGs in OG')
    plt.ylabel('Number of OGs')
    plt.title('Distribution of OGs across number of pOGs')
    plt.legend()
    plt.savefig(f'out/{k}clique/OGnum-pOGnum_type_dist2-1.png')
    plt.xlim((-1, 17))  # Adjust axis to truncate outliers
    plt.savefig(f'out/{k}clique/OGnum-pOGnum_type_dist2-2.png')
    plt.close()

    plt.pie([sum(OGtype.values()) for OGtype in OGtypes], labels=[f'Type {i}' for i in range(len(OGtypes))])
    plt.title('OGs by type')
    plt.savefig(f'out/{k}clique/type_pie.png')
    plt.close()

    print()
    print(f'{k}-CLIQUE')
    for i, OGtype in enumerate(OGtypes):
        print(f'Type {i}:', sum(OGtype.values()))


# Load seq metadata
gnid2ppids = {}
ppid2gnid = {}
ppid2spid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, repr = line.split()
        ppid2spid[ppid] = spid
        ppid2gnid[ppid] = gnid
        if repr == 'True':
            try:
                gnid2ppids[gnid].append(ppid)
            except KeyError:
                gnid2ppids[gnid] = [ppid]

# Load pgraph
pgraph = {}
with open('../hits2pgraph/out/pgraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        pgraph[node] = bitscores

# Load OGs
OGs = {}
with open('../clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        OGs[OGid] = gnids

ks = list(range(4, 7))
pOGs_ks = {k: [] for k in ks}
OGtypes_ks = {k: [{} for _ in range(5)] for k in ks}
for OGid, OG in OGs.items():
    # Make subpgraph
    subpgraph = {}
    ppids = set([ppid for gnid in OG for ppid in gnid2ppids[gnid]])
    for ppid in ppids:
        # Collect hits by SPID
        spids = {}
        for adj_node, adj_bitscore in pgraph.get(ppid, {}).items():  # In case PPID has no hits
            spid = ppid2spid[adj_node]
            try:
                spids[spid][adj_node] = adj_bitscore
            except KeyError:
                spids[spid] = {adj_node: adj_bitscore}

        # Find max hits for each SPID
        s = set()
        for spid, adjs in spids.items():
            max_bitscore = max(adjs.values())
            s.update([adj_node for adj_node, adj_bitscore in adjs.items() if adj_bitscore == max_bitscore])
        subpgraph[ppid] = s

    # Filter subpgraph by reciprocity
    for node, adjs in subpgraph.items():
        del_elems = []
        for adj in adjs:
            if not (adj in subpgraph and node in subpgraph[adj]):
                del_elems.append(adj)
        for del_elem in del_elems:
            adjs.remove(del_elem)

    # Create graph
    G = nx.Graph()
    for node, adjs in subpgraph.items():
        G.add_node(node)
        for adj in adjs:
            G.add_edge(node, adj)

    # Handle cliques
    try:
        signal.signal(signal.SIGALRM, percolate.clique_handler)
        if OGid in ['096d', '09e0', '09dd', '09e1', '09de']:
            raise percolate.CliqueError
        signal.alarm(30)
        cliques = list(find_cliques(G))
        signal.alarm(0)
    except percolate.CliqueError:
        print(f'CliqueError: {OGid}')
        for k in ks:
            subpOGs = set()
            core = k_core(G, k)
            for component in connected_components(core):
                subpOGs.add(frozenset([frozenset(edge) for edge in core.edges(component)]))
            pOGs_ks[k].append(subpOGs)
            classify_OG(OGtypes_ks[k], subpOGs)
        continue  # Continue to next OG

    # Handle percolation
    for k in ks:
        try:
            signal.signal(signal.SIGALRM, percolate.percolate_handler)
            signal.alarm(60)
            subpOGs = list(percolate.k_clique_communities_progressive(G, k, cliques))
            signal.alarm(0)
        except percolate.PercolateError:
            print(f'PercolateError: ({k}, {OGid})')
            subpOGs = set()
            core = k_core(G, k)
            for component in connected_components(core):
                subpOGs.add(frozenset([frozenset(edge) for edge in core.edges(component)]))
        pOGs_ks[k].append(subpOGs)
        classify_OG(OGtypes_ks[k], subpOGs)

for k in ks:
    save_results(pOGs_ks[k], OGtypes_ks[k], k)

"""
CliqueError: 096d
CliqueError: 09dd
CliqueError: 09de
CliqueError: 09e0
CliqueError: 09e1
CliqueError: 3020

4-CLIQUE
Type 0: 25
Type 1: 12085
Type 2: 150
Type 3: 1170
Type 4: 1103

5-CLIQUE
Type 0: 100
Type 1: 11966
Type 2: 322
Type 3: 1085
Type 4: 1060

6-CLIQUE
Type 0: 1015
Type 1: 10272
Type 2: 1047
Type 3: 1059
Type 4: 1140

DEPENDENCIES
../../ortho_search/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../clique4+_gcommunity/out/clique4+_gcommunity2.py
    ../clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt
../hits2pgraph/hits2pgraph.py
    ../hits2pgraph/out/pgraph2.tsv
"""