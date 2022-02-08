"""Convert HSPs to hits."""

import multiprocessing as mp
import numpy as np
import os
from itertools import groupby


def parse_file(qspid, sspid):
    # Open files and process lines
    with open(f'../blast2hsps/out/hsps/{qspid}/{sspid}') as file:
        file.readline()  # Skip header

        hits = []
        for _, group in groupby(file, line2key):
            hsps = []
            for line in group:
                hsp = {column: f(field) for (column, f), field in zip(hsp_columns.items(), line.split())}
                if hsp['bitscore'] >= 50:
                    hsps.append(hsp)
            if hsps:
                hits.append(hsps2hit(hsps))

    # Make output directories
    if not os.path.exists(f'out/{qspid}/'):
        os.makedirs(f'out/{qspid}/')  # Recursive folder creation

    # Write lines to output
    with open(f'out/{qspid}/{sspid}', 'w') as file:
        file.write('\t'.join(hit_columns) + '\n')
        for hit in hits:
            file.write('\t'.join([str(hit[column]) for column in hit_columns]) + '\n')


def hsps2hit(hsps):
    # Calculate values from disjoint HSPs only
    output_hsps = [hsp for hsp in hsps if hsp['disjoint']]
    hit = {**{key: output_hsps[0][key] for key in ['qppid', 'qgnid', 'qspid', 'sppid', 'sgnid', 'sspid', 'qlen', 'slen']},
           'hspnum': len(output_hsps),
           'nqa': sum([hsp['qend'] - hsp['qstart'] + 1 for hsp in output_hsps]),
           'nsa': sum([hsp['send'] - hsp['sstart'] + 1 for hsp in output_hsps]),
           'bitscore': sum([hsp['bitscore'] for hsp in output_hsps])}

    # Calculate values from disjoint and compatible HSPs
    nondisjoint_hsps = [hsp for hsp in hsps if not hsp['disjoint']]
    for hsp in nondisjoint_hsps:
        if is_compatible(hsp, output_hsps, 'query') and is_compatible(hsp, output_hsps, 'subject'):
            output_hsps.append(hsp)
    hit['chspnum'] = len(output_hsps)
    qcov = np.zeros((1, output_hsps[0]['qlen']), dtype=bool)
    scov = np.zeros((1, output_hsps[0]['slen']), dtype=bool)
    for hsp in output_hsps:
        qcov[0, hsp['qstart']-1:hsp['qend']] = True
        scov[0, hsp['sstart']-1:hsp['send']] = True
    hit['cnqa'] = qcov.sum()
    hit['cnsa'] = scov.sum()

    return hit


def line2key(line):
    fields = line.rstrip('\n').split('\t')
    return fields[0], fields[3]


def is_compatible(hsp, hsp_list, key_type):
    """Return if hsp is compatible with all HSPs in hsp_list.

    For each test_hsp in hsp_list, overlap of hsp and test_hsp must not exceed
    50% of the length of either. If so, returns True, else returns False.
    """
    if key_type == 'query':
        start_key = 'qstart'
        end_key = 'qend'
    elif key_type == 'subject':
        start_key = 'sstart'
        end_key = 'send'
    else:
        raise ValueError('key_type is not query or subject')

    for test_hsp in hsp_list:
        if not (hsp[start_key] > test_hsp[end_key] or test_hsp[start_key] > hsp[end_key]):
            start = max(hsp[start_key], test_hsp[start_key])
            end = min(hsp[end_key], test_hsp[end_key])
            length = end-start
            if length / (hsp[end_key]-hsp[start_key]) >= 0.5 or length / (test_hsp[end_key]-test_hsp[start_key]) >= 0.5:
                return False
    return True


hsp_columns = {'qppid': str, 'qgnid': str, 'qspid': str,
               'sppid': str, 'sgnid': str, 'sspid': str,
               'length': int, 'nident': int, 'gaps': int,
               'qlen': int, 'qstart': int, 'qend': int,
               'slen': int, 'sstart': int, 'send': int,
               'evalue': float, 'bitscore': float,
               'index_hsp': lambda x: x == 'True',
               'disjoint': lambda x: x == 'True'}
hit_columns = ['qppid', 'qgnid', 'qspid',
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
        pool.starmap(parse_file, spids)

"""
DEPENDENCIES
../blast2hsps/blast2hsps.py
    ../blast2hsps/out/hsps/*/*.tsv
"""