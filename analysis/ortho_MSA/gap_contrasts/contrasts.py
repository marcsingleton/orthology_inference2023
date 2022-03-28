"""Calculate gap contrasts between alignments."""

import os
import re

import numpy as np
import pandas as pd
import skbio
from src.utils import read_fasta


def get_contrasts(node):
    child1, child2 = node.children
    if child1.is_tip():
        contrasts1, val1, bl1 = [], child1.value, child1.length
    else:
        contrasts1, val1, bl1 = get_contrasts(child1)
    if child2.is_tip():
        contrasts2, val2, bl2 = [], child2.value, child2.length
    else:
        contrasts2, val2, bl2 = get_contrasts(child2)

    bl_sum = bl1 + bl2
    value = (val1 * bl2 + val2 * bl1) / bl_sum
    branch_length = node.length + bl1 * bl2 / bl_sum
    contrasts = contrasts1 + contrasts2
    contrasts.append((val1 - val2) / bl_sum)

    return contrasts, value, branch_length


# Load tree
tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
spids = {tip.name for tip in tree_template.tips() if tip.name != 'sleb'}

# Load representative OGs
OG_filter = pd.read_table('../OG_filter/out/OG_filter.tsv')

# Calculate contrasts
totals, rows = [], []
for record in OG_filter.itertuples():
    if record.ppdidnum == record.gnidnum:
        msa = read_fasta(f'../align_fastas1/out/{record.OGid}.afa')
    else:
        msa = read_fasta(f'../align_fastas2/out/{record.OGid}.afa')
    msa = {re.search(r'spid=([a-z]+)', header).group(1): seq for header, seq in msa}

    tree = tree_template.deepcopy().shear(msa.keys())
    for tip in tree.tips():
        gap_vector = np.asarray([1 if sym == '-' else 0 for sym in msa[tip.name]])
        tip.value = gap_vector
    tree.length = 0  # Set root length to 0 for convenience

    contrasts, _, _ = get_contrasts(tree)
    gap_matrix = np.asarray([[0 if sym == '-' else 1 for sym in seq] for seq in msa.values()])
    len1 = len(msa['dmel'])  # Total length of alignment
    len2 = (gap_matrix / len(msa)).sum()  # Adjusted length of alignment
    totals.append([record.OGid, str(len(msa)), str(len1), str(len2), str(np.abs(contrasts).sum())])
    if len(msa) == len(spids):
        row_sums = list(np.abs(contrasts).sum(axis=1))
        rows.append([record.OGid, str(len1), str(len2)] + [str(row_sum) for row_sum in row_sums])

# Write contrasts to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/total_sums.tsv', 'w') as file:
    header = '\t'.join(['OGid', 'gnidnum', 'len1', 'len2', 'total']) + '\n'
    file.write(header)
    for total in totals:
        file.write('\t'.join(total) + '\n')
with open('out/row_sums.tsv', 'w') as file:
    header = '\t'.join(['OGid', 'len1', 'len2'] + [f'row{i}' for i in range(len(spids)-1)]) + '\n'
    file.write(header)
    for row in rows:
        file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.afa
../align_fastas2/align_fastas2.py
    ../align_fastas2/out/*.afa
../OG_filter/OG_filter.py
    ../OG_filter/out/OG_filter.tsv
"""