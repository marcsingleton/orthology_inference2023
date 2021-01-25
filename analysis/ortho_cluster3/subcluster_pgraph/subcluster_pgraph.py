"""Extract clusters from OGs."""

import matplotlib.pyplot as plt
import os
from itertools import combinations
from src.ortho_cluster.triDFS import cluster

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
OGs = []
with open('../clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        _, _, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        OGs.append(gnids)

pOGs = []
OGtypes = [{} for _ in range(5)]
for OG in OGs:
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

    # Cluster by triangle criterion
    subpOGs = cluster(subpgraph)
    pOGs.append(subpOGs)

    # Classify OGs
    subnOGs = [set([node for edge in subpOG for node in edge]) for subpOG in subpOGs]
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

# Make plots output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write pOGs to file
j = 0
with open('out/pclusters.txt', 'w') as file:
    for i, subpOGs in enumerate(pOGs):
        OGid = hex(i)[2:].zfill(4)
        for pOG in subpOGs:
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
plt.savefig(f'out/OGnum-pOGnum_type_dist1-1.png')
plt.xlim((-1, 17))  # Adjust axis to truncate outliers
plt.savefig(f'out/OGnum-pOGnum_type_dist1-2.png')
plt.close()

plt.bar(OGtypes[3].keys(), OGtypes[3].values(), label='Type 3', color='C3')
plt.bar(OGtypes[4].keys(), OGtypes[4].values(), bottom=[OGtypes[3].get(key, 0) for key in OGtypes[4]], label='Type 4', color='C4')
plt.xlabel('Number of pOGs in OG')
plt.ylabel('Number of OGs')
plt.title('Distribution of OGs across number of pOGs')
plt.legend()
plt.savefig(f'out/OGnum-pOGnum_type_dist2-1.png')
plt.xlim((-1, 17))  # Adjust axis to truncate outliers
plt.savefig(f'out/OGnum-pOGnum_type_dist2-2.png')
plt.close()

plt.pie([sum(OGtype.values()) for OGtype in OGtypes], labels=[f'Type {i}' for i in range(len(OGtypes))])
plt.title('OGs by type')
plt.savefig(f'out/type_pie.png')
plt.close()

for i, OGtype in enumerate(OGtypes):
    print(f'Type {i}:', sum(OGtype.values()))

"""
OUTPUT
Type 0: 2
Type 1: 11987
Type 2: 68
Type 3: 1414
Type 4: 1062

DEPENDENCIES
../../../src/ortho_cluster/triDFS.py
../../ortho_search/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../clique4+_gcommunity/out/clique4+_gcommunity2.py
    ../clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt
../hits2pgraph/hits2pgraph.py
    ../hits2pgraph/out/pgraph2.tsv
"""