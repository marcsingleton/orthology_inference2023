"""Extract metadata from OGs including numbers of edges, genes, species, and sequences."""

import os
import pandas as pd

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2meta[ppid] = (gnid, spid, sqid)

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

# Load gOGs
OGid2gOGid = {}
with open('../../ortho_cluster3/connect_OGgraph/out/OGconnect.txt') as file:
    for line in file:
        gOGid, OGids = line.rstrip().split(':')
        for OGid in OGids.split(','):
            OGid2gOGid[OGid] = gOGid

# Load OGs
rows = []
with open('../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        gnids = set([ppid2meta[ppid][0] for ppid in ppids])
        spids = set([ppid2meta[ppid][1] for ppid in ppids])
        sqids = set([ppid2meta[ppid][2] for ppid in ppids])
        bitscore = 0
        for edge in edges.split('\t'):
            node1, node2 = edge.split(',')
            bitscore += pgraph[node1][node2] + pgraph[node2][node1]

        d = {'CCid': CCid, 'OGid': OGid, 'gOGid': OGid2gOGid[OGid],
             'bitscore': round(bitscore, 1), 'edgenum': len(edges.split('\t')),
             'ppidnum': len(ppids), 'sqidnum': len(sqids),
             'gnidnum': len(gnids), 'spidnum': len(spids)}
        rows.append(d)
OGs = pd.DataFrame(rows)

# Find representatives
idxs = []
for _, group in OGs.groupby('gOGid'):
    idx = max(group.itertuples(), key=lambda x: (x.spidnum, -x.gnidnum, x.bitscore)).Index  # Most species, least genes, highest bitscore
    idxs.append(idx)
OGs['repr'] = False
OGs.loc[idxs, 'repr'] = True
repr_OGs = OGs[OGs['repr']]

# Print counts
num = 31

ppnum = OGs['ppidnum'] == num
sqnum = OGs['sqidnum'] == num
gnnum = OGs['gnidnum'] == num
spnum = OGs['spidnum'] == num

print('ALL OGS')
print('Total OGs:', len(OGs))
print(f'OGs with {num} genes:', len(OGs[gnnum]))
print(f'OGs with {num} genes and species:', len(OGs[gnnum & spnum]))
print(f'OGs with {num} genes, species, and sequences:', len(OGs[gnnum & spnum & ppnum]))
print(f'OGs with {num} genes, species, and unique sequences:', len(OGs[gnnum & spnum & sqnum]))

ppnum = repr_OGs['ppidnum'] == num
sqnum = repr_OGs['sqidnum'] == num
gnnum = repr_OGs['gnidnum'] == num
spnum = repr_OGs['spidnum'] == num

print()
print('REPRESENTATIVE OGS')
print('Total OGs:', len(repr_OGs))
print(f'Representative OGs with {num} genes:', len(repr_OGs[gnnum]))
print(f'Representative OGs with {num} genes and species:', len(repr_OGs[gnnum & spnum]))
print(f'Representative OGs with {num} genes, species, and sequences:', len(repr_OGs[gnnum & spnum & ppnum]))
print(f'Representative OGs with {num} genes, species, and unique sequences:', len(repr_OGs[gnnum & spnum & sqnum]))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

OGs.to_csv('out/OG_meta.tsv', sep='\t', index=False)

"""
OUTPUT 
ALL OGS
Total OGs: 23283
OGs with 31 genes: 7703
OGs with 31 genes and species: 7487
OGs with 31 genes, species, and sequences: 2047
OGs with 31 genes, species, and unique sequences: 5091

REPRESENTATIVE OGS
Total OGs: 15099
Representative OGs with 31 genes: 7484
Representative OGs with 31 genes and species: 7290
Representative OGs with 31 genes, species, and sequences: 1942
Representative OGs with 31 genes, species, and unique sequences: 4931

DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity2.py
    ../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../../ortho_cluster3/connect_OGgraph/connect_OGgraph.py
    ../../ortho_cluster3/connect_OGgraph/out/OGconnect.txt
../../ortho_cluster3/hits2pgraph/hits2pgraph.py
    ../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv
"""