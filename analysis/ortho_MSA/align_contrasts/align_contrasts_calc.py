"""Calculate gap contrasts between alignments."""

import os
import re

import numpy as np
import pandas as pd
import skbio


def get_contrasts(node):
    child1, child2 = node.children

    if child1.is_tip() and child2.is_tip():
        contrasts1, val1, bl1 = [], child1.value, child1.length
        contrasts2, val2, bl2 = [], child2.value, child2.length
    elif child1.is_tip():
        contrasts1, val1, bl1 = [], child1.value, child1.length
        contrasts2, val2, bl2 = get_contrasts(child2)
    elif child2.is_tip():
        contrasts1, val1, bl1 = get_contrasts(child1)
        contrasts2, val2, bl2 = [], child2.value, child2.length
    else:
        contrasts1, val1, bl1 = get_contrasts(child1)
        contrasts2, val2, bl2 = get_contrasts(child2)

    bl_sum = bl1 + bl2
    value = (val1 * bl2 + val2 * bl1) / bl_sum
    branch_length = node.length + bl1 * bl2 / bl_sum
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


# Load tree
tree_template = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree_template = tree_template.shear([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])

# Load pOGs and extract representatives for each OG
pOG_meta = pd.read_table('../pOG_meta/out/pOG_meta.tsv')

idxs = []
for _, group in pOG_meta.groupby(['OGid', 'overlap']):
    idx = max(group.itertuples(), key=lambda x: (x.spidnum, -x.gnidnum, x.bitscore)).Index  # Most species, least genes, highest bitscore
    idxs.append(idx)
df1 = pOG_meta.iloc[idxs]
df2 = df1[(df1['gnidnum'] == 26) & (df1['spidnum'] == 26)]

# Calculate contrasts
if not os.path.exists('out/'):
    os.mkdir('out/')

rows = []
for record in df2.itertuples():
    if record.ppidnum == record.gnidnum:
        MSA = load_alignment(f'../align_fastas1/out/{record.pOGid}.mfa')
    else:
        MSA = load_alignment(f'../align_fastas2-2/out/{record.pOGid}.mfa')

    tree = tree_template.deepcopy()
    for tip in tree.tips():
        gap_vector = np.asarray([1 if sym == '-' else 0 for sym in MSA[tip.name]])
        tip.value = gap_vector
    tree.length = 0  # Set root length to 0 for convenience

    contrasts, _, _ = get_contrasts(tree)
    row_sums = list(np.abs(contrasts).sum(axis=1))
    gap_matrix = np.asarray([[0 if sym == '-' else 1 for sym in seq] for seq in MSA.values()])
    len1 = len(MSA['dmel'])  # Total length of alignment
    len2 = (gap_matrix / 26).sum()  # Adjusted length of alignment
    rows.append([record.pOGid, str(len1), str(len2)] + [str(row_sum) for row_sum in row_sums])

with open('out/row_sums.tsv', 'w') as file:
    header = '\t'.join(['pOGid', 'len1', 'len2'] + [f'c{i}' for i in range(25)]) + '\n'
    file.writelines(header)
    for row in rows:
        file.writelines('\t'.join(row) + '\n')

"""
DEPENDENCIES
../ortho_tree/consensus_tree/consensus_tree.py
    ../ortho_tree/consensus_tree/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../pOG_meta/pOG_meta.py
    ../pOG_meta/out/pOG_meta.tsv
"""