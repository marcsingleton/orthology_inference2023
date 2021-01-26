"""Extract metadata from pOGs including numbers of edges, genes, species, and sequences."""

import os
from itertools import groupby

import pandas as pd

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
with open('../../ortho_cluster3/clique4+_pcommunity/out/5clique/pclusters.txt') as file:
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
            ds[pOGid] = {'OGid': OGid, 'pOGid': pOGid, 'overlap': False, 'bitscore': round(bitscore, 1),
                         'edgenum': len(edges.split('\t')), 'ppidnum': len(ppids), 'gnidnum': len(gnids), 'spidnum': len(spids)}
        for pOGid in pOGid2gnids:  # Calculate overlap
            gnids1 = pOGid2gnids[pOGid]
            gnids2 = set().union(*[value for key, value in pOGid2gnids.items() if key != pOGid])
            ds[pOGid]['overlap'] = bool(gnids1 & gnids2)
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
Total OGs: 20685
OGs with 26 genes: 7938
OGs with 26 genes and species: 7748
OGs with 26 genes, species, and sequences: 5663

DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity.py
    ../../ortho_cluster3/clique4+_pcommunity/out/5clique/pclusters.txt
../../ortho_cluster3/hits2pgraph/hits2pgraph.py
    ../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv
"""