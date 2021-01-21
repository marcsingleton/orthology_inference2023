"""Plot distribution of bitscores associated with edges in clustered ggraph."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from numpy import linspace

# Load seq metadata
gnid2spid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        _, gnid, spid, _ = line.split()
        gnid2spid[gnid] = spid

# Load ggraph
ggraph = {}
with open('../hits2ggraph/out/ggraph1.tsv') as file:
    for line in file:
        node1, adjs = line.rstrip('\n').split('\t')
        bitscores = [adj.split(':') for adj in adjs.split(',')]
        ggraph[node1] = {node2: int(bitscore) for (node2, bitscore) in bitscores}

# Extract OG hits
rows = []
with open('../subcluster_ggraph/out/ggraph1/gclusters.txt') as file:
    for line in file:
        _, _, edges = line.rstrip().split(':')
        if len(edges.split('\t')) > 3:  # Discard edges in minimal OGs
            for edge in edges.split('\t'):
                node1, node2 = edge.split(',')
                row1 = {'qgnid': node1, 'sgnid': node2,
                        'qspid': gnid2spid[node1], 'sspid': gnid2spid[node2],
                        'bitscore': ggraph[node1][node2]}
                row2 = {'qgnid': node2, 'sgnid': node1,
                        'qspid': gnid2spid[node2], 'sspid': gnid2spid[node1],
                        'bitscore': ggraph[node2][node1]}
                rows.append(row1)
                rows.append(row2)
hits = pd.DataFrame(rows)

# Make output directory
if not os.path.exists('out/ggraph1/'):
    os.makedirs('out/ggraph1/')  # Recursive folder creation

bins = linspace(0, hits['bitscore'].max(), 300, endpoint=True)

# 1 Query plots
if not os.path.exists('out/ggraph1/qhists/'):
    os.mkdir('out/ggraph1/qhists/')
g = hits.groupby('qspid')

means = g['bitscore'].mean().sort_index()
plt.bar(means.index, means.values, width=0.75)
plt.xticks(rotation=60)
plt.xlabel('Species')
plt.ylabel('Mean bitscore as query')
plt.savefig('out/ggraph1/bar_bitscoremean-qspecies.png')
plt.close()

for qspid, group in g:
    plt.hist(group['bitscore'], bins=bins)
    plt.xlabel('Bitscore')
    plt.ylabel('Number of hits')
    plt.savefig(f'out/ggraph1/qhists/{qspid}.png')
    plt.close()

# 2 Subject plots
if not os.path.exists('out/ggraph1/shists/'):
    os.mkdir('out/ggraph1/shists/')
g = hits.groupby('sspid')

means = g['bitscore'].mean().sort_index()
plt.bar(means.index, means.values, width=0.75)
plt.xticks(rotation=60)
plt.xlabel('Species')
plt.ylabel('Mean bitscore as subject')
plt.savefig('out/ggraph1/bar_bitscoremean-sspecies.png')
plt.close()

for sspid, group in g:
    plt.hist(group['bitscore'], bins=bins)
    plt.xlabel('Bitscore')
    plt.ylabel('Number of hits')
    plt.savefig(f'out/ggraph1/shists/{sspid}.png')
    plt.close()

# 3 Query-subject plots
if not os.path.exists('out/ggraph1/qshists/'):
    os.mkdir('out/ggraph1/qshists/')
g = hits.groupby(['qspid', 'sspid'])

means = g['bitscore'].mean()
plt.hist(means, bins=50)
plt.xlabel('Mean bitscore')
plt.ylabel('Number of query-subject pairs')
plt.savefig('out/ggraph1/hist_qsnum-bitscoremean.png')
plt.close()

for (qspid, sspid), group in g:
    plt.hist(group['bitscore'], bins=bins)
    plt.xlabel('Bitscore')
    plt.ylabel('Number of hits')
    plt.savefig(f'out/ggraph1/qshists/{qspid}-{sspid}.png')
    plt.close()

"""
DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../hits2ggraph/out/hits2ggraph1.py
    ../hits2ggraph/out/ggraph1.tsv
../subcluster_ggraph/subcluster_ggraph1.py
    ../subcluster_ggraph/out/ggraph1/gclusters.txt
"""