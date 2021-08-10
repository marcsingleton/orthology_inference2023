"""Convert hits to an undirected pgraph."""

import multiprocessing as mp
import os
import pandas as pd
from itertools import permutations


def load_hit(qspid, sspid):
    df = pd.read_csv(f'../../ortho_search/hsps2hits/out/{qspid}/{sspid}.tsv', sep='\t',
                     usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
    r = pd.read_csv(f'../../ortho_search/hits2reciprocal/out/{qspid}/{sspid}.tsv', sep='\t',
                    usecols=['reciprocal1'], memory_map=True)

    return df[r['reciprocal1']]


dtypes = {'qppid': 'string', 'sppid': 'string',
          'bitscore': float}
num_processes = 2

if __name__ == '__main__':
    # Parse genomes
    spids = []
    with open('../config/genomes.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            spids.append(line.split()[0])

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

    # Make output directory
    if not os.path.exists('out/'):
        os.mkdir('out/')

    # Write to file
    with open('out/pgraph1.tsv', 'w') as file:
        for qgnid, edges in graph.items():
            file.write(qgnid + '\t' + ','.join([sgnid + ':' + str(bitscore) for sgnid, bitscore in edges]) + '\n')

"""
DEPENDENCIES
../../ortho_search/hsps2hits/hsps2hits.py
    ../../ortho_search/hsps2hits/out/*/*.tsv
../../ortho_search/hits2reciprocal/hits2reciprocal.py
    ../../ortho_search/hits2reciprocal/out/*/*.tsv
../config/genomes.tsv
"""