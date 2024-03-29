"""Convert hits to an undirected graph."""

import multiprocessing as mp
import os
import pandas as pd
from itertools import permutations


def load_hit(qspid, sspid):
    df = pd.read_table(f'../../ortho_search/hsps2hits/out/{qspid}/{sspid}.tsv',
                       usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
    r = pd.read_table(f'../../ortho_search/hits2reciprocal/out/{qspid}/{sspid}.tsv',
                      usecols=['reciprocal'], memory_map=True)

    return df[r['reciprocal']]


dtypes = {'qppid': 'string', 'sppid': 'string',
          'bitscore': float}
num_processes = 2

if __name__ == '__main__':
    # Load genomes
    spids = []
    with open('../config/genomes.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            spids.append(fields['spid'])

    # Load data
    with mp.Pool(processes=num_processes) as pool:
        hits = pd.concat(pool.starmap(load_hit, permutations(spids, 2)))

    graph = {}
    for row in hits.itertuples():
        qppid, sppid, bitscore = row.qppid, row.sppid, row.bitscore
        try:
            graph[qppid].append((sppid, float(bitscore)))
        except KeyError:
            graph[qppid] = [(sppid, float(bitscore))]

    # Write to file
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open('out/hit_graph.tsv', 'w') as file:
        for qppid, edges in graph.items():
            file.write(qppid + '\t' + ','.join([sppid + ':' + str(bitscore) for sppid, bitscore in edges]) + '\n')

"""
DEPENDENCIES
../../ortho_search/hsps2hits/hsps2hits.py
    ../../ortho_search/hsps2hits/out/*/*.tsv
../../ortho_search/hits2reciprocal/hits2reciprocal.py
    ../../ortho_search/hits2reciprocal/out/*/*.tsv
../config/genomes.tsv
"""