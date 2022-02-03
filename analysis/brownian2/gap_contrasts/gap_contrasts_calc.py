"""Calculate gap contrasts between alignments."""

import os
import re

import numpy as np
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


ppid_regex = r'ppid=([A-Za-z0-9_]+)'
spid_regex = r'spid=([a-z]+)'

# Load tree
tree_template = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
spids = set([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])

# Load regions
regions = []
with open('../aucpred_filter/out/regions_30.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.split()
        regions.append((OGid, int(start), int(stop), set(ppids.split(','))))

# Calculate contrasts
if not os.path.exists('out/'):
    os.mkdir('out/')

totals = []
rows = []
for OGid, start, stop, ppids in regions:
    msa = {}
    for header, seq in read_fasta(f'../insertion_trim/out/{OGid}.mfa'):
        ppid = re.search(ppid_regex, header).group(1)
        spid = re.search(spid_regex, header).group(1)
        if ppid in ppids:
            msa[spid] = seq[start:stop]

    tree = tree_template.deepcopy().shear(msa.keys())
    for tip in tree.tips():
        gap_vector = np.asarray([1 if sym == '-' else 0 for sym in msa[tip.name]])
        tip.value = gap_vector
    tree.length = 0  # Set root length to 0 for convenience

    contrasts, _, _ = get_contrasts(tree)
    gap_matrix = np.asarray([[0 if sym == '-' else 1 for sym in seq] for seq in msa.values()])
    len1 = len(msa['dmel'])  # Total length of alignment
    len2 = (gap_matrix / len(msa)).sum()  # Adjusted length of alignment
    totals.append([OGid, str(start), str(stop), str(len(msa)), str(len1), str(len2), str(np.abs(contrasts).sum())])
    if len(msa) == len(spids):
        row_sums = list(np.abs(contrasts).sum(axis=1))
        rows.append([OGid, str(start), str(stop), str(len1), str(len2)] + [str(row_sum) for row_sum in row_sums])

with open('out/total_sums.tsv', 'w') as file:
    header = '\t'.join(['OGid', 'start', 'stop', 'gnidnum', 'len1', 'len2', 'total']) + '\n'
    file.writelines(header)
    for total in totals:
        file.writelines('\t'.join(total) + '\n')
with open('out/row_sums.tsv', 'w') as file:
    header = '\t'.join(['OGid', 'start', 'stop', 'len1', 'len2'] + [f'row{i}' for i in range(len(spids)-1)]) + '\n'
    file.writelines(header)
    for row in rows:
        file.writelines('\t'.join(row) + '\n')

"""
DEPENDENCIES
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
../aucpred_filter/aucpred_filter.py
    ../aucpred_filter/out/regions_30.tsv
"""