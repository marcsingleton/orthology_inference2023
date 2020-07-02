"""Convert OG graph into tabular format."""

import json
import os
import pandas as pd


def extract_edges(qgnid, sgnid, ggraph, gnid2spid, OGid):
    qspid, sspid = gnid2spid[qgnid], gnid2spid[sgnid]
    qppids = ggraph[qgnid][sgnid]

    rows = []
    for qppid, sppids in qppids.items():
        for sppid, params in sppids.items():
            row = {'qppid': qppid, 'qgnid': qgnid, 'qspid': qspid,
                   'sppid': sppid, 'sgnid': sgnid, 'sspid': sspid,
                   'OGid': OGid, **params}
            rows.append(row)
    return rows


# Load gn metadata
gnid2spid = {}
with open('../ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        _, gnid, spid = line.split()
        gnid2spid[gnid] = spid

# Load best hits graph
with open('../make_xreciprocal/out/ggraph.json') as file:
    ggraph = json.load(file)

# Load OGs
OGs = {}
with open('../subcluster_xgraph/out/ggraph/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        OGs[OGid] = [tuple(edge.split(',')) for edge in edges.split('\t')]

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Extract edge data from graph
rows = []
for OGid, OG in OGs.items():
    for node1, node2 in OG:
        rows.extend(extract_edges(node1, node2, ggraph, gnid2spid, OGid))
        rows.extend(extract_edges(node2, node1, ggraph, gnid2spid, OGid))
df = pd.DataFrame(rows)
df.to_csv('out/ggraph.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../make_xreciprocal/make_greciprocal.py
    ../make_xreciprocal/out/ggraph.json
../ppid2meta/ppid2meta.py
    ../ppid2meta/out/ppid2meta.tsv
../subcluster_xgraph/subcluster_ggraph.py
    ../subcluster_xgraph/out/ggraph/gclusters.txt
"""