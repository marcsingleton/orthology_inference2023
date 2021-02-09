"""Extract metadata from pOGs including numbers of edges, genes, species, and sequences."""

import os
from itertools import combinations, groupby

import pandas as pd
from src.ortho_cluster.DFS import connect

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, _ = line.split()
        ppid2meta[ppid] = (gnid, spid)

# Load pgraph
pgraph = {}
with open('../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        pgraph[node] = bitscores

# Load pOGs
rows = []
with open('../../ortho_cluster3/subcluster_pgraph/out/pclusters.txt') as file:
    for OGid, pOGs in groupby(file, lambda line: line.rstrip().split(':')[0]):
        pOGid2gnids = {}  # Store gnids separately to test for overlaps
        ds = {}  # Store rows by pOGid to modify repr after all GNIDs are collected
        for line in pOGs:  # Calculate metadata
            _, pOGid, edges = line.rstrip().split(':')
            ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
            gnids = set([ppid2meta[ppid][0] for ppid in ppids])
            spids = set([ppid2meta[ppid][1] for ppid in ppids])
            bitscore = 0
            for edge in edges.split('\t'):
                node1, node2 = edge.split(',')
                bitscore += pgraph[node1][node2] + pgraph[node2][node1]

            pOGid2gnids[pOGid] = gnids
            ds[pOGid] = {'OGid': OGid, 'pOGid': pOGid, 'pCCid': None, 'bitscore': round(bitscore, 1),
                         'edgenum': len(edges.split('\t')), 'ppidnum': len(ppids), 'gnidnum': len(gnids), 'spidnum': len(spids)}

        # Calculate overlaps
        graph = {pOGid: set() for pOGid in pOGid2gnids}
        for pOGid1, pOGid2 in combinations(pOGid2gnids, 2):
            gnids1 = pOGid2gnids[pOGid1]
            gnids2 = pOGid2gnids[pOGid2]
            if gnids1 & gnids2:
                graph[pOGid1].add(pOGid2)
                graph[pOGid2].add(pOGid1)
        pCCs = connect(graph)
        for i, pCC in enumerate(pCCs):
            pCCid = hex(i)[2:].zfill(4)
            for pOGid in pCC:
                ds[pOGid]['pCCid'] = pCCid
        rows.extend(ds.values())
pOGs = pd.DataFrame(rows)

# Print counts
num = 26
ppnum = pOGs['ppidnum'] == num
gnnum = pOGs['gnidnum'] == num
spnum = pOGs['spidnum'] == num

print('Total pOGs:', len(pOGs))
print(f'pOGs with {num} genes:', len(pOGs[gnnum]))
print(f'pOGs with {num} genes and species:', len(pOGs[gnnum & spnum]))
print(f'pOGs with {num} genes, species, and sequences:', len(pOGs[gnnum & spnum & ppnum]))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

pOGs.to_csv('out/pOG_meta.tsv', sep='\t', index=False)

"""
OUTPUT 
Total pOGs: 21632
pOGs with 26 genes: 8156
pOGs with 26 genes and species: 7970
pOGs with 26 genes, species, and sequences: 5256

DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_cluster3/subcluster_pgraph/subcluster_pgraph.py
    ../../ortho_cluster3/subcluster_pgraph/out/pclusters.txt
../../ortho_cluster3/hits2pgraph/hits2pgraph.py
    ../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv
"""