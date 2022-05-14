"""Convert HSPs to hits."""

import multiprocessing as mp
import os
from itertools import groupby, permutations

import numpy as np


def parse_file(qspid, sspid):
    # Open files and process lines
    with open(f'../blast2hsps/out/hsps/{qspid}/{sspid}.tsv') as file:
        file.readline()  # Skip header

        hits = []
        for _, group in groupby(file, line2key):
            hsps = []
            for line in group:
                hsp = {column: f(field) for (column, f), field in zip(hsp_columns.items(), line.rstrip('\n').split('\t'))}
                if hsp['compatible']:  # Use compatible HSPs only (which pass the E-value cutoff)
                    hsps.append(hsp)
            if hsps:
                hits.append(hsps2hit(hsps))

    # Write to file
    if not os.path.exists(f'out/{qspid}/'):
        os.makedirs(f'out/{qspid}/')

    with open(f'out/{qspid}/{sspid}.tsv', 'w') as file:
        file.write('\t'.join(hit_columns) + '\n')
        for hit in hits:
            file.write('\t'.join([str(hit[column]) for column in hit_columns]) + '\n')


def hsps2hit(hsps):
    # Calculate values from all HSPs
    hit = {key: hsps[0][key] for key in ['qppid', 'qgnid', 'sppid', 'sgnid', 'qlen', 'slen']}
    hit['chspnum'] = len(hsps)
    qcov = np.zeros((1, hsps[0]['qlen']), dtype=bool)
    scov = np.zeros((1, hsps[0]['slen']), dtype=bool)
    for hsp in hsps:
        qcov[0, hsp['qstart']-1:hsp['qend']] = True
        scov[0, hsp['sstart']-1:hsp['send']] = True
    hit['cnqa'] = qcov.sum()
    hit['cnsa'] = scov.sum()

    # Calculate values from disjoint HSPs only
    disjoint_hsps = [hsp for hsp in hsps if hsp['disjoint']]
    hit['hspnum'] = len(disjoint_hsps)
    hit['nqa'] = sum([hsp['qend'] - hsp['qstart'] + 1 for hsp in disjoint_hsps])
    hit['nsa'] = sum([hsp['send'] - hsp['sstart'] + 1 for hsp in disjoint_hsps])
    hit['bitscore'] = sum([hsp['bitscore'] for hsp in disjoint_hsps])

    return hit


def line2key(line):
    fields = line.rstrip('\n').split('\t')
    return fields[0], fields[2]


hsp_columns = {'qppid': str, 'qgnid': str,
               'sppid': str, 'sgnid': str,
               'length': int, 'nident': int, 'gaps': int,
               'qlen': int, 'qstart': int, 'qend': int,
               'slen': int, 'sstart': int, 'send': int,
               'evalue': float, 'bitscore': float,
               'index_hsp': lambda x: x == 'True',
               'disjoint': lambda x: x == 'True',
               'compatible': lambda x: x == 'True'}
hit_columns = ['qppid', 'qgnid',
               'sppid', 'sgnid',
               'hspnum', 'chspnum',
               'qlen', 'nqa', 'cnqa',
               'slen', 'nsa', 'cnsa',
               'bitscore']
num_processes = 2

# Load genomes
spids = []
with open('../config/genomes.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        spids.append(fields['spid'])

# Parse HSPs
if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(parse_file, permutations(spids, 2))

"""
DEPENDENCIES
../config/genomes.tsv
../blast2hsps/blast2hsps.py
    ../blast2hsps/out/hsps/*/*.tsv
"""