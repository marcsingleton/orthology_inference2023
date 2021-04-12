"""Trim alignments, removing spurious sequences and regions."""

import json
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import skbio
from src.ortho_MSA.trim import trim_conserved, trim_insertions


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


def trim(OGid):
    # 0 Load MSA
    try:
        msa1 = load_msa(f'../align_fastas1/out/{OGid}.mfa')
    except FileNotFoundError:
        msa1 = load_msa(f'../align_fastas2-2/out/{OGid}.mfa')

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
        msa2.append((seq.metadata['description'], ''.join(syms)))

    # 4 Remove gap only columns
    slices, idx = [], None
    for j in range(len(msa2[0][1])):
        for i in range(len(msa2)):
            sym = msa2[i][1][j]
            if sym != '-':
                if idx is None:  # Store position only if new slice is not started
                    idx = j
                break
        else:
            if idx is not None:
                slices.append(slice(idx, j))
                idx = None
    if idx is not None:  # Add final slice to end
        slices.append(slice(idx, len(msa2[0][1])))

    # 5 Write to file
    with open(f'out/{OGid}.mfa', 'w') as file:
        for header, seq1 in msa2:
            seq2 = []
            for s in slices:
                seq2.extend(seq1[s])
            seq2 = ''.join(seq2)

            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)]) + '\n'
            file.write(header)
            file.write(seqstring)


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