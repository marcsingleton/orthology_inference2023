"""Trim alignments, removing spurious sequences and regions."""

import json
import os

import pandas as pd
import skbio
from src.ortho_MSA.trim import trim_msa

# Load parameters
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

for row in OG_filter.itertuples():
    # Load MSA
    if row.sqidnum == row.gnidnum:
        msa1 = skbio.read(f'../align_fastas1/out/{row.OGid}.mfa',
                          format='fasta', into=skbio.TabularMSA, constructor=skbio.Protein)
    else:
        msa1 = skbio.read(f'../align_fastas2-2/out/{row.OGid}.mfa',
                          format='fasta', into=skbio.TabularMSA, constructor=skbio.Protein)

    # Trim MSA and remove gap only columns
    msa2 = trim_msa(msa1,
                    tp['con_frac'], tp['con_window'], tp['con_minlen'], tp['con_rate'], tp['con_minsig'],
                    tp['gap_num'], tp['gap_rate'], tp['gap_minsig'],
                    tp['gp_sigma'], tp['gd_window'], tp['indel1_rate'], tp['indel2_rate'],
                    tp['weights'], tp['threshold'],
                    matrix)
    slices, i = [], None
    for j, col in enumerate(msa2.iter_positions()):
        for sym in col:
            if sym != '-':
                if i is None:  # Store position only if new slice is not started
                    i = j
                break
        else:
            slices.append(slice(i, j))
            i = None
    if i is not None:  # Add final slice to end
        slices.append(slice(i, msa2.shape[1]))
    msa3 = msa2.loc[:, slices]

    skbio.write(msa3, format='fasta', into=f'out/{row.OGid}.mfa', max_width=80)

"""
DEPENDENCIES
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../OG_filter/OG_filter.py
    ../OG_filter/out/OG_filter.tsv
"""