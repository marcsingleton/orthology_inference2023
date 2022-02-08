"""Compute reciprocity for each hit."""

import os


def check_edge(qppid, sppid, graph):
    # Check both directions since not all hits are in graph
    try:
        r1 = qppid in graph[sppid]
        r2 = sppid in graph[qppid]
        r = r1 and r2
    except KeyError:
        r = False
    return r


columns = {'qppid': str, 'qgnid': str, 'qspid': str,
           'sppid': str, 'sgnid': str, 'sspid': str,
           'hspnum': int, 'chspnum': int,
           'qlen': int, 'nqa': int, 'cnqa': int,
           'slen': int, 'nsa': int, 'cnsa': int,
           'bitscore': float}

# Load graphs
graph1 = {}
with open('../hits2pgraph/out/pgraph1.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph1[node] = set([adj.split(':')[0] for adj in adjs.split(',')])
graph2 = {}
with open('../hits2pgraph/out/pgraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph2[node] = set([adj.split(':')[0] for adj in adjs.split(',')])

# Check reciprocity
for qspid in os.listdir('../hsps2hits/out/'):
    for sspid in os.listdir(f'../hsps2hits/out/{qspid}/'):
        rows = []
        with open(f'../hsps2hits/out/{qspid}/{sspid}') as file:
            file.readline()  # Skip header
            for line in file:
                hit = {column: f(field) for (column, f), field in zip(columns.items(), line.split())}
                qppid, sppid = hit['qppid'], hit['sppid']
                qlen, cnqa = hit['qlen'], hit['cnqa']

                r1 = check_edge(qppid, sppid, graph1)
                r2 = cnqa / qlen >= 0.5 and check_edge(qppid, sppid, graph2)
                rows.append((qppid, sppid, str(r1), str(r2)))

        # Make output directory
        if not os.path.exists(f'out/{qspid}/'):
            os.makedirs(f'out/{qspid}/')  # Recursive folder creation

        # Write to file
        with open(f'out/{qspid}/{sspid}', 'w') as file:
            file.write('qppid\tsppid\treciprocal1\treciprocal2\n')
            for row in rows:
                file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../hsps2hits/hsps2hits.py
    ../hsps2hits/out/*/*.tsv
../hits2pgraph/hits2pgraph.py
    ../hits2pgraph/out/pgraph.tsv
"""