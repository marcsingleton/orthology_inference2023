"""Extract hits which are contained in OGs."""

import multiprocessing as mp
import os
import pandas as pd


def get_ids(qgnid, sgnid, edge2ids):
    try:
        return edge2ids[(qgnid, sgnid)]
    except KeyError:
        try:
            return edge2ids[(sgnid, qgnid)]
        except KeyError:
            return 'NA', 'NA'


def parse_hit(qspid, sspid):
    df = pd.read_csv(f'../hsps2hits/out/{qspid}/{sspid}', sep='\t',
                     usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
    r = pd.read_csv(f'../hits2reciprocal/out/{qspid}/{sspid}', sep='\t',
                    usecols=['reciprocal1', 'reciprocal2'], memory_map=True)
    dfr = df.join(r)

    orows = []
    gnids = set()
    for irow in dfr.itertuples():
        qgnid, sgnid = irow.qgnid, irow.sgnid
        if (qgnid, sgnid) in gnids:
            continue

        CCid1, OGid1 = get_ids(qgnid, sgnid, edge2ids1)
        CCid2, OGid2 = get_ids(qgnid, sgnid, edge2ids2)

        orows.append((qgnid, sgnid, CCid1, OGid1, CCid2, OGid2))
        gnids.add((qgnid, sgnid))

    # Make output directory
    if not os.path.exists(f'out/{qspid}/'):
        os.makedirs(f'out/{qspid}/')  # Recursive folder creation

    # Write to file
    with open(f'out/{qspid}/{sspid}', 'w') as file:
        file.write('qgnid\tsgnid\tCCid1\tOGid1\tCCid2\tOGid2\n')
        for orow in orows:
            file.write('\t'.join(orow) + '\n')


dtypes = {'qgnid': 'string', 'sgnid': 'string'}
num_processes = 2

if __name__ == '__main__':
    # Load OG edges
    edge2ids1 = {}
    with open('../subcluster_ggraph/out/ggraph1/gclusters.txt') as file:
        for line in file:
            CCid, OGid, edges = line.rstrip().split(':')
            edge2ids1.update({tuple(edge.split(',')): (CCid, OGid) for edge in edges.split('\t')})
    edge2ids2 = {}
    with open('../subcluster_ggraph/out/ggraph2/gclusters.txt') as file:
        for line in file:
            _, _, edges = line.rstrip().split(':')
            CCid, OGid, edges = line.rstrip().split(':')
            edge2ids2.update({tuple(edge.split(',')): (CCid, OGid) for edge in edges.split('\t')})

    # Iterate through hits
    with mp.Pool(processes=num_processes) as pool:
        tsvs = [(qspid, sspid) for qspid in os.listdir('../hsps2hits/out/')
                for sspid in os.listdir(f'../hsps2hits/out/{qspid}/')]
        pool.starmap(parse_hit, tsvs)

"""
DEPENDENCIES
../hits2reciprocal/hits2reciprocal.py
    ../hits2reciprocal/out/*/*.tsv
../hsps2hits/hsps2hits.py
    ../hsps2hits/out/*/*.tsv
../subcluster_ggraph/subcluster_ggraph1.py
    ../subcluster_ggraph/out/ggraph1/gclusters.txt
../subcluster_ggraph/subcluster_ggraph2.py
    ../subcluster_ggraph/out/ggraph2/gclusters.txt
"""