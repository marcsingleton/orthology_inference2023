"""Filter OGs to representatives by uniqueness criteria and species number."""

import os

import matplotlib.pyplot as plt
import pandas as pd


def spid_filter(spids):
    conditions = [({'dnov', 'dvir'}, 1),
                  ({'dmoj', 'dnav'}, 1),
                  ({'dinn', 'dgri', 'dhyd'}, 2),
                  ({'dgua', 'dsob'}, 1),
                  ({'dbip', 'dana'}, 1),
                  ({'dser', 'dkik'}, 1),
                  ({'dele', 'dfik'}, 1),
                  ({'dtak', 'dbia'}, 1),
                  ({'dsuz', 'dspu'}, 1),
                  ({'dsan', 'dyak'}, 1),
                  ({'dmel'}, 1),
                  ({'dmau', 'dsim', 'dsec'}, 1)]
    return all([len(spids & group) >= num for group, num in conditions])


# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2meta[ppid] = (gnid, spid, sqid)

# Load graph
graph = {}
with open('../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        graph[node] = bitscores

# Load gOGs
OGid2gOGid = {}
with open('../../ortho_cluster3/connect_OGgraph/out/OGconnect.txt') as file:
    for line in file:
        gOGid, OGids = line.rstrip().split(':')
        for OGid in OGids.split(','):
            OGid2gOGid[OGid] = gOGid

# Load OGs
gOGid2OGs = {}
with open('../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = {node for edge in edges.split('\t') for node in edge.split(',')}
        gnids = {ppid2meta[ppid][0] for ppid in ppids}
        spids = {ppid2meta[ppid][1] for ppid in ppids}
        sqids = {ppid2meta[ppid][2] for ppid in ppids}
        bitscore = 0
        for edge in edges.split('\t'):
            node1, node2 = edge.split(',')
            bitscore += graph[node1][node2] + graph[node2][node1]

        if len(gnids) == len(spids) >= 20 and spid_filter(spids):
            gOGid = OGid2gOGid[OGid]
            record = {'CCid': CCid, 'OGid': OGid, 'gOGid': gOGid,
                      'bitscore': round(bitscore, 1), 'edgenum': len(edges.split('\t')),
                      'ppidnum': len(ppids), 'sqidnum': len(sqids),
                      'gnidnum': len(gnids), 'spidnum': len(spids)}
            try:
                gOGid2OGs[gOGid].append(record)
            except KeyError:
                gOGid2OGs[gOGid] = [record]

rows = []
for OGs in gOGid2OGs.values():
    OG = max(OGs, key=lambda x: (x['spidnum'], x['bitscore']))  # Most species, highest bitscore
    rows.append(OG)
OGs = pd.DataFrame(rows)

if not os.path.exists('out/'):
    os.mkdir('out/')

print('Total filtered OGs:', len(OGs))
counts = OGs['gnidnum'].value_counts()
plt.bar(counts.index, counts.values)
plt.xlabel('Number of genes in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/bar_OGnum-gnnum.png')

OGs.to_csv('out/OG_filter.tsv', sep='\t', index=False)

"""
OUTPUT
Total filtered OGs: 8551

DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity2.py
    ../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../../ortho_cluster3/connect_OGgraph/connect_OGgraph.py
    ../../ortho_cluster3/connect_OGgraph/out/OGconnect.txt
../../ortho_cluster3/hits2pgraph/hits2pgraph2.py
    ../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv
"""