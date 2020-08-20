"""Convert hits to an undirected ggraph."""

import multiprocessing as mp
import os
import pandas as pd


def load_hit(qspid, sspid):
    df = pd.read_csv(f'../hsps2hits/out/{qspid}/{sspid}', sep='\t',
                     usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
    r = pd.read_csv(f'../hits2reciprocal/out/{qspid}/{sspid}', sep='\t',
                    usecols=['reciprocal1', 'reciprocal2'], memory_map=True)
    dfr = df.join(r)

    return dfr[dfr['reciprocal1'] | dfr['reciprocal2']]


def add_edge(qgnid, sgnid, r, ggraph):
    if r:
        try:
            ggraph[qgnid].add(sgnid)
        except KeyError:
            ggraph[qgnid] = set([sgnid])


dtypes = {'qppid': 'string', 'qgnid': 'string',
          'sppid': 'string', 'sgnid': 'string'}
num_processes = 2

if __name__ == '__main__':
    # Load data
    with mp.Pool(processes=num_processes) as pool:
        tsvs = [(qspid, sspid) for qspid in os.listdir('../hsps2hits/out/')
                for sspid in os.listdir(f'../hsps2hits/out/{qspid}/')]
        hits = pd.concat(pool.starmap(load_hit, tsvs))

    ggraph1 = {}
    ggraph2 = {}
    for row in hits.itertuples():
        qgnid, sgnid = row.qgnid, row.sgnid
        r1, r2 = row.reciprocal1, row.reciprocal2

        add_edge(qgnid, sgnid, r1, ggraph1)
        add_edge(qgnid, sgnid, r2, ggraph2)

    # Make output directory
    if not os.path.exists('out/'):
        os.mkdir('out/')

    # Write to file
    with open('out/ggraph1.tsv', 'w') as file:
        for qgnid, sgnids in ggraph1.items():
            file.write(qgnid + '\t' + ','.join(sgnids) + '\n')
    with open('out/ggraph2.tsv', 'w') as file:
        for qgnid, sgnids in ggraph2.items():
            file.write(qgnid + '\t' + ','.join(sgnids) + '\n')

"""
DEPENDENCIES
../hsps2hits/hsps2hits.py
    ../hsps2hits/out/*/*.tsv
../hits2reciprocal/hits2reciprocal.py
    ../hits2reciprocal/out/*/*.tsv
"""