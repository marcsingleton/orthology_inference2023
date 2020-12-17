"""Calculate gap contrasts between alignments."""

import os
import re
from copy import deepcopy

import Bio.Phylo as Phylo
import numpy as np
import pandas as pd


def get_contrasts(node):
    leaf1, leaf2 = node.clades

    if leaf1.is_terminal() and leaf2.is_terminal():
        contrasts1, val1, bl1 = [], leaf1.value, leaf1.branch_length
        contrasts2, val2, bl2 = [], leaf2.value, leaf2.branch_length
    elif leaf1.is_terminal():
        contrasts1, val1, bl1 = [], leaf1.value, leaf1.branch_length
        contrasts2, val2, bl2 = get_contrasts(leaf2)
    elif leaf2.is_terminal():
        contrasts1, val1, bl1 = get_contrasts(leaf1)
        contrasts2, val2, bl2 = [], leaf2.value, leaf2.branch_length
    else:
        contrasts1, val1, bl1 = get_contrasts(leaf1)
        contrasts2, val2, bl2 = get_contrasts(leaf2)

    bl_sum = bl1 + bl2
    value = (val1 * bl2 + val2 * bl1) / bl_sum
    branch_length = node.branch_length + bl1 * bl2 / bl_sum
    contrasts = contrasts1 + contrasts2
    contrasts.append((val1 - val2) / bl_sum)

    return contrasts, value, branch_length


def load_alignment(path):
    MSA = {}
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                spid = re.search(r'spid=([a-z]+)', line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            MSA[spid] = seq
    return MSA


tree_template = Phylo.read('../consensus_tree/out/100red_ni.txt', 'newick')
tree_template.root_at_midpoint()  # Removes multifurcation at dwil

OGid2meta = pd.read_table('../OGid2meta/out/OGid2meta.tsv').drop(['CCid', 'edgenum'], axis=1)
df = OGid2meta[(OGid2meta['gnidnum'] == 25) & (OGid2meta['spidnum'] == 25)]

if not os.path.exists('out/'):
    os.makedirs('out/contrasts')

rows = []
for record in df.itertuples():
    if record.sqidnum == 25:
        MSA = load_alignment(f'../align_fastas1/out/{record.OGid}.mfa')
    else:
        MSA = load_alignment(f'../align_fastas2/out/{record.OGid}.mfa')

    tree = deepcopy(tree_template)
    for terminal in tree.get_terminals():
        gap_vector = np.asarray([1 if sym == '-' else 0 for sym in MSA[terminal.name]])
        terminal.value = gap_vector

    contrasts, _, _ = get_contrasts(tree.root)
    contrasts = np.asarray(contrasts)  # Convert to 2D array
    with open(f'out/contrasts/{record.OGid}.npy', 'wb') as file:
        np.save(file, contrasts)

    row_sums = list(np.abs(contrasts).sum(axis=1))
    gap_matrix = np.asarray([[0 if sym == '-' else 1 for sym in seq] for seq in MSA.values()])
    len1 = len(MSA['dmel'])  # Total length of alignment
    len2 = (gap_matrix / 25).sum()  # Adjusted length of alignment
    rows.append([record.OGid, str(len1), str(len2)] + [str(row_sum) for row_sum in row_sums])

with open('out/row_sums.tsv', 'w') as file:
    header = '\t'.join(['OGid', 'len1', 'len2'] + [f'c{i}' for i in range(24)]) + '\n'
    file.writelines(header)
    for row in rows:
        file.writelines('\t'.join(row) + '\n')

"""
DEPENDENCIES
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2/align_fastas2.py
    ../align_fastas2/out/*.mfa
../consensus_tree/consensus_tree.py
    ../consensus_tree/out/100red_ni.txt
../OGid2meta/OGid2meta.py
    ../OGid2meta/out/OGid2meta.tsv
"""