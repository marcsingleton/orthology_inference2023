"""Trim alignments, removing spurious sequences and regions."""

import json
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import skbio
from src.ortho_MSA.trim import trim_conserved, trim_insertions
from src.utils import read_fasta


def trim(OGid):
    # 0 Load MSA
    try:
        msa1 = read_fasta(f'../align_fastas1/out/{OGid}.mfa')
    except FileNotFoundError:
        msa1 = read_fasta(f'../align_fastas2-2/out/{OGid}.mfa')

    # 1 Calculate shared variables
    gaps_array = np.full((len(msa1), len(msa1[0][1])), False)
    for i, (_, seq) in enumerate(msa1):
        for j, sym in enumerate(seq):
            if sym == '-':
                gaps_array[i, j] = True
    scores = gaps_array.sum(axis=0)
    msa1 = skbio.TabularMSA([skbio.Protein(seq, metadata={'description': header}) for header, seq in msa1])

    # 2 Get trims (segments and columns)
    syms_list1 = trim_conserved(msa1, scores, gaps_array,
                                tp['con_frac'], tp['con_window'], tp['con_minlen'], tp['con_rate'], tp['con_minsig'])
    syms_list2, trims = trim_insertions(msa1, scores, gaps_array,
                                        tp['gap_num'], tp['gap_rate'], tp['gap_minsig'],
                                        tp['nongap_frac'], tp['nongap_minlen'],
                                        tp['gp_sigma'], tp['gd_window'], tp['indel1_rate'], tp['indel2_rate'],
                                        tp['weights'], tp['threshold'],
                                        matrix)

    # 3 Combine trims (segments and columns) to yield final alignment
    msa2 = []
    for seq, syms1, syms2 in zip(msa1, syms_list1, syms_list2):
        syms = ['-' if sym1 != sym2 else sym1 for sym1, sym2 in zip(syms1, syms2)]  # Will only differ if one is converted to gap
        msa2.append((seq.metadata['description'], syms))

    # 4 Restore gap only columns
    gaps_array = np.full((len(msa2), len(msa2[0][1])), False)
    for i, (_, seq) in enumerate(msa2):
        for j, sym in enumerate(seq):
            if sym == '-':
                gaps_array[i, j] = True
    scores = gaps_array.sum(axis=0)

    rf = ['x' for _ in range(len(msa2[0][1]))]  # Metadata for marking consensus columns in profile HMM
    for region, in ndimage.find_objects(ndimage.label(scores == len(msa2))[0]):
        rf[region] = (region.stop - region.start) * ['.']
        for i in range(len(msa2)):
            syms = msa2[i][1]
            syms[region] = list(str(msa1[i, region]))

    # 5 Write to file
    msa2 = skbio.TabularMSA([skbio.Protein(''.join(syms), metadata={'description': header}) for header, syms in msa2],
                            positional_metadata={'RF': rf})
    msa2.write(f'out/{OGid}.sto', 'stockholm')


# Load parameters
num_processes = int(os.environ['SLURM_NTASKS'])

with open('trim_params.json') as file:
    tp = json.load(file)  # Trim parameters
with open('../config/trim_params.json') as file:
    tp.update(json.load(file))
matrix = {}
with open('../config/BLOSUM62.txt') as file:
    file.readline()  # Skip header
    syms = file.readline().split()
    for i, line in enumerate(file):
        for j, value in enumerate(line.split()[1:]):
            matrix[(syms[i], syms[j])] = int(value)

OG_filter = pd.read_table('../OG_filter/out/OG_filter.tsv')

if not os.path.exists('out/'):
    os.mkdir('out/')

if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        pool.map(trim, OG_filter['OGid'])

"""
DEPENDENCIES
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../OG_filter/OG_filter.py
    ../OG_filter/out/OG_filter.tsv
"""