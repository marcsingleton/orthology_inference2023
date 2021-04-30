"""Convert HSPs to hits."""

import multiprocessing as mp
import numpy as np
import os
from itertools import groupby


def parse_hsp(qspid, sspid):
    # Open files and process lines
    with open(f'../blast2hsps/out/hsps/{qspid}/{sspid}') as file:
        file.readline()

        hits = []
        for _, group in groupby(file, line2key):
            ds = []
            for line in group:
                d = {column: f(field) for (column, f), field in zip(icolumns.items(), line.split())}
                if d['bitscore'] >= 50:
                    ds.append(d)
            if ds:
                hits.append(hsps2hit(ds))

    # Make output directories
    if not os.path.exists(f'out/{qspid}/'):
        os.makedirs(f'out/{qspid}/')  # Recursive folder creation

    # Write lines to output
    with open(f'out/{qspid}/{sspid}', 'w') as file:
        file.write('\t'.join(ocolumns) + '\n')
        for hit in hits:
            file.write('\t'.join([str(hit[column]) for column in ocolumns]) + '\n')


def hsps2hit(hsps):
    # Extract parameters from disjoint HSPs
    ohsps = [hsp for hsp in hsps if hsp['disjoint']]
    row = {**{key: ohsps[0][key] for key in ['qppid', 'qgnid', 'qspid', 'sppid', 'sgnid', 'sspid', 'qlen', 'slen']},
           'hspnum': len(ohsps),
           'nqa': sum([ohsp['qend'] - ohsp['qstart'] + 1 for ohsp in ohsps]),
           'nsa': sum([ohsp['send'] - ohsp['sstart'] + 1 for ohsp in ohsps]),
           'bitscore': sum([ohsp['bitscore'] for ohsp in ohsps])}

    # Extract parameters from non-compatible HSPs
    ihsps = [hsp for hsp in hsps if not hsp['disjoint']]
    for ihsp in ihsps:
        if not has_overlap(ihsp, ohsps, True) and not has_overlap(ihsp, ohsps, False):
            ohsps.append(ihsp)
    row['chspnum'] = len(ohsps)
    qcov = np.zeros((1, ohsps[0]['qlen']), dtype=bool)
    scov = np.zeros((1, ohsps[0]['slen']), dtype=bool)
    for ohsp in ohsps:
        qcov[0, ohsp['qstart']-1:ohsp['qend']] = True
        scov[0, ohsp['sstart']-1:ohsp['send']] = True
    row['cnqa'] = qcov.sum()
    row['cnsa'] = scov.sum()

    return row


def line2key(line):
    fields = line.split()
    return fields[0], fields[3]


def has_overlap(ihsp, ohsps, query=True):
    if query:
        start_key = 'qstart'
        end_key = 'qend'
    else:
        start_key = 'sstart'
        end_key = 'send'

    for ohsp in ohsps:
        if not (ihsp[start_key] > ohsp[end_key] or ohsp[start_key] > ihsp[end_key]):
            start = max(ihsp[start_key], ohsp[start_key])
            end = min(ihsp[end_key], ohsp[end_key])
            if (end-start) / (ihsp[end_key]-ihsp[start_key]) >= 0.5 or (end-start) / (ohsp[end_key]-ohsp[start_key]) >= 0.5:
                return True
    return False


icolumns = {'qppid': str, 'qgnid': str, 'qspid': str,
            'sppid': str, 'sgnid': str, 'sspid': str,
            'length': int, 'nident': int, 'gaps': int,
            'qlen': int, 'qstart': int, 'qend': int,
            'slen': int, 'sstart': int, 'send': int,
            'evalue': float, 'bitscore': float,
            'index_hsp': lambda x: True if x == 'True' else False,
            'disjoint': lambda x: True if x == 'True' else False}
ocolumns = ['qppid', 'qgnid', 'qspid',
            'sppid', 'sgnid', 'sspid',
            'hspnum', 'chspnum',
            'qlen', 'nqa', 'cnqa',
            'slen', 'nsa', 'cnsa',
            'bitscore']
num_processes = 2

if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        spids = [(qspid, sspid) for qspid in os.listdir('../blast2hsps/out/hsps/')
                 for sspid in os.listdir(f'../blast2hsps/out/hsps/{qspid}')]
        pool.starmap(parse_hsp, spids)

"""
DEPENDENCIES
../blast2hsps/blast2hsps.py
    ../blast2hsps/out/hsps/*/*.tsv
"""