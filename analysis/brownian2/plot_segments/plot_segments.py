"""Plot ground-truth segment labels."""

import os

import matplotlib.pyplot as plt
import numpy as np
from src.brownian2.trim import trim_terminals
from src.draw import plot_msa_lines


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


OGid2labels = {}
with open('../config/segments.tsv') as file:
    file.readline()
    for line in file:
        fields = line.split()
        OGid, start, stop, state = fields
        if OGid in OGid2labels:
            OGid2labels[OGid][state].append((int(start), int(stop)))
        else:
            d = {'0': [], '1A': [], '1B': [], '2': [], '3': []}
            d[state].append((int(start), int(stop)))
            OGid2labels[OGid] = d

if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid, labels in OGid2labels.items():
    msa = trim_terminals(load_msa(f'../../ortho_MSA/realign_hmmer/out/{OGid}.mfa'))

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

    plot_msa_lines([seq[1].upper() for seq in msa], [lines['1A'], lines['2'], lines['3'], lines['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}.png', bbox_inches='tight')
    plt.close()

"""
DEPENDENCIES
../config/segments.tsv
"""