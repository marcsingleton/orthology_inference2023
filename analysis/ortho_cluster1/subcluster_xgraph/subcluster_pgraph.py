"""Extract clusters from connected components of pgraph."""

import matplotlib.pyplot as plt
import os
from itertools import combinations
from src.ortho_cluster.graphs import get_triangle_clusters

# Load pgraph
pgraph = {}
with open('../hsps2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        pgraph[node] = adjs.split(',')

# Make reciprocal
rpgraph = {}
for qppid, sppids in pgraph.items():
    for sppid in sppids:
        try:
            r = qppid in pgraph[sppid]
        except KeyError:
            r = False

        if r:
            try:
                rpgraph[qppid].add(sppid)
            except KeyError:
                rpgraph[qppid] = {sppid}

# Load connected components
CCs = []
with open('../connect_xgraph/out/pconnect.txt') as file:
    for line in file:
        _, nodes = line.rstrip().split(':')
        CCs.append(set(nodes.split(',')))

OGs = []
CCtypes = [{} for _ in range(5)]
for CC in CCs:
    subpgraph = {node: rpgraph[node] for node in CC}

    # Cluster by triangle criterion
    subOGs = get_triangle_clusters(subpgraph)
    OGs.append(subOGs)

    # Classify CCs
    subnOGs = [{node for edge in subOG for node in edge} for subOG in subOGs]
    if len(subnOGs) == 0:
        CCtypes[0][len(subnOGs)] = CCtypes[0].get(len(subnOGs), 0) + 1  # Component has no OGs
    elif len(subnOGs) == 1:
        if len(subnOGs[0]) == len(CC):
            CCtypes[1][len(subnOGs)] = CCtypes[1].get(len(subnOGs), 0) + 1  # Component and OG are equal
        else:
            CCtypes[2][len(subnOGs)] = CCtypes[2].get(len(subnOGs), 0) + 1  # Component has single OG which is a subset of the component
    elif any([set.intersection(nOG1, nOG2) for nOG1, nOG2 in combinations(subnOGs, 2)]):
        CCtypes[3][len(subnOGs)] = CCtypes[4].get(len(subnOGs), 0) + 1  # Component has multiple non-disjoint OGs
    else:
        CCtypes[4][len(subnOGs)] = CCtypes[3].get(len(subnOGs), 0) + 1  # Component has multiple pairwise disjoint OGs

# Make plots output directory
if not os.path.exists('out/pgraph/'):
    os.makedirs('out/pgraph/')  # Recursive folder creation

# Write OGs to file
j = 0
with open('out/pgraph/pclusters.txt', 'w') as file:
    for i, subOGs in enumerate(OGs):
        CCid = hex(i)[2:].zfill(4)
        for OG in subOGs:
            OGid = hex(j)[2:].zfill(4)
            file.write(CCid + ':' + OGid + ':' + '\t'.join([f'{node1},{node2}' for node1, node2 in OG]) + '\n')
            j += 1

# Plots
plt.bar(CCtypes[0].keys(), CCtypes[0].values(), label='Type 0')
plt.bar(CCtypes[1].keys(), CCtypes[1].values(), label='Type 1')
plt.bar(CCtypes[2].keys(), CCtypes[2].values(), bottom=CCtypes[1][1], label='Type 2')
plt.bar(CCtypes[3].keys(), CCtypes[3].values(), label='Type 3')
plt.bar(CCtypes[4].keys(), CCtypes[4].values(), bottom=[CCtypes[3].get(key, 0) for key in CCtypes[4]], label='Type 4')
plt.xlabel('Number of OGs in connected component')
plt.ylabel('Number of connected components')
plt.title('Distribution of connected components across number of OGs')
plt.legend()
plt.savefig('out/pgraph/connectnum-OGnum_type_dist1.png')
plt.close()

plt.bar(CCtypes[3].keys(), CCtypes[3].values(), label='Type 3', color='C3')
plt.bar(CCtypes[4].keys(), CCtypes[4].values(), bottom=[CCtypes[3].get(key, 0) for key in CCtypes[4]], label='Type 4', color='C4')
plt.xlabel('Number of OGs in connected component')
plt.ylabel('Number of connected components')
plt.title('Distribution of connected components across number of OGs')
plt.legend()
plt.savefig('out/pgraph/connectnum-OGnum_type_dist2.png')
plt.close()

plt.pie([sum(CCtype.values()) for CCtype in CCtypes], labels=[f'Type {i}' for i in range(len(CCtypes))])
plt.title('Connected Components by type')
plt.savefig('out/pgraph/type_pie.png')
plt.close()

for i, CCtype in enumerate(CCtypes):
    print(f'Type {i}:', sum(CCtype.values()))

"""
OUTPUT
Type 0: 17602
Type 1: 11236
Type 2: 1328
Type 3: 243
Type 4: 237

DEPENDENCIES
../connect_xgraph/connect_pgraph.py
    ../connect_xgraph/out/pconnect.txt
../hsps2pgraph/hsps2pgraph.py
    ../hsps2pgraph/out/pgraph.tsv
"""