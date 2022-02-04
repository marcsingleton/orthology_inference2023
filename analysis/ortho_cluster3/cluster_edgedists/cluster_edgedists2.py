"""Plot distribution of bitscores associated with edges in clustered pgraph."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from numpy import linspace

# Load seq metadata
ppid2spid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, _, spid, _ = line.split()
        ppid2spid[ppid] = spid

# Load pgraph
graph = {}
with open('../hits2pgraph/out/pgraph2.tsv') as file:
    for line in file:
        node1, adjs = line.rstrip('\n').split('\t')
        bitscores = [adj.split(':') for adj in adjs.split(',')]
        graph[node1] = {node2: float(bitscore) for (node2, bitscore) in bitscores}

# Extract OG hits
rows = []
with open('../clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        _, _, edges = line.rstrip().split(':')
        if len(edges.split('\t')) > 4:  # Discard edges in minimal OGs
            for edge in edges.split('\t'):
                node1, node2 = edge.split(',')
                row1 = {'qppid': node1, 'sppid': node2,
                        'qspid': ppid2spid[node1], 'sspid': ppid2spid[node2],
                        'bitscore': graph[node1][node2]}
                row2 = {'qppid': node2, 'sppid': node1,
                        'qspid': ppid2spid[node2], 'sspid': ppid2spid[node1],
                        'bitscore': graph[node2][node1]}
                rows.append(row1)
                rows.append(row2)
hits = pd.DataFrame(rows)

# Make output directory
if not os.path.exists('out/pgraph2/'):
    os.makedirs('out/pgraph2/')  # Recursive folder creation

bins = linspace(0, hits['bitscore'].max(), 300, endpoint=True)

# 1 Query plots
if not os.path.exists('out/pgraph2/qhists/'):
    os.mkdir('out/pgraph2/qhists/')
g = hits.groupby('qspid')

means = g['bitscore'].mean().sort_index()
plt.bar(means.index, means.values, width=0.75)
plt.xticks(rotation=60, fontsize=8)
plt.xlabel('Species')
plt.ylabel('Mean bitscore as query')
plt.subplots_adjust(bottom=0.15)
plt.savefig('out/pgraph2/bar_bitscoremean-qspecies.png')
plt.close()

for qspid, group in g:
    plt.hist(group['bitscore'], bins=bins)
    plt.xlabel('Bitscore')
    plt.ylabel('Number of hits')
    plt.savefig(f'out/pgraph2/qhists/{qspid}.png')
    plt.close()

# 2 Subject plots
if not os.path.exists('out/pgraph2/shists/'):
    os.mkdir('out/pgraph2/shists/')
g = hits.groupby('sspid')

means = g['bitscore'].mean().sort_index()
plt.bar(means.index, means.values, width=0.75)
plt.xticks(rotation=60, fontsize=8)
plt.xlabel('Species')
plt.ylabel('Mean bitscore as subject')
plt.subplots_adjust(bottom=0.15)
plt.savefig('out/pgraph2/bar_bitscoremean-sspecies.png')
plt.close()

for sspid, group in g:
    plt.hist(group['bitscore'], bins=bins)
    plt.xlabel('Bitscore')
    plt.ylabel('Number of hits')
    plt.savefig(f'out/pgraph2/shists/{sspid}.png')
    plt.close()

# 3 Query-subject plots
if not os.path.exists('out/pgraph2/qshists/'):
    os.mkdir('out/pgraph2/qshists/')
g = hits.groupby(['qspid', 'sspid'])

means = g['bitscore'].mean()
plt.hist(means, bins=50)
plt.xlabel('Mean bitscore')
plt.ylabel('Number of query-subject pairs')
plt.savefig('out/pgraph2/hist_qsnum-bitscoremean.png')
plt.close()

for (qspid, sspid), group in g:
    plt.hist(group['bitscore'], bins=bins)
    plt.xlabel('Bitscore')
    plt.ylabel('Number of hits')
    plt.savefig(f'out/pgraph2/qshists/{qspid}-{sspid}.png')
    plt.close()

"""
DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../clique4+_pcommunity/clique4+_pcommunity2.py
    ../clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../hits2pgraph/hits2pgraph2.py
    ../hits2pgraph/out/pgraph2.tsv
"""