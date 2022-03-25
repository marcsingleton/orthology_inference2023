"""Plot ground-truth segment labels."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import skbio
from src.brownian2.trim import trim_terminals
from src.draw import plot_msa_lines
from src.utils import read_fasta

OGid2labels = {}
with open('../config/segments.tsv') as file:
    file.readline()
    for line in file:
        fields = line.rstrip('\n').split('\t')
        OGid, start, stop, state = fields
        if OGid in OGid2labels:
            OGid2labels[OGid][state].append((int(start), int(stop)))
        else:
            labels = {'0': [], '1A': [], '1B': [], '2': [], '3': []}
            labels[state].append((int(start), int(stop)))
            OGid2labels[OGid] = labels

tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}

if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid, labels in OGid2labels.items():
    msa = trim_terminals(read_fasta(f'../../ortho_MSA/realign_hmmer/out/{OGid}.afa'))
    msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]

    if labels['0'] and labels['0'][0][0] == 0:
        offset = labels['0'][0][1]
    else:
        offset = 0

    lines = {}
    for state in ['1A', '1B', '2', '3']:
        line = np.zeros(len(msa[0][1]))
        for start, stop in labels[state]:
            line[slice(start-offset, stop-offset)] = 1
        lines[state] = line

    msa = [seq.upper() for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only
    plot_msa_lines(msa, [lines['1A'], lines['2'], lines['3'], lines['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}.png', bbox_inches='tight')
    plt.close()

"""
DEPENDENCIES
../../ortho_MSA/realign_hmmer/realign_hmmer.py
    ../../ortho_MSA/realign_hmmer/out/*.afa
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../config/segments.tsv
"""