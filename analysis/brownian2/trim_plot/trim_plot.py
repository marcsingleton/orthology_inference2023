"""Extract state 1 segments to yield trimmed alignments."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio
from src.draw import plot_msa_lines
from src.brownian2.trim import trim_terminals, get_slices


def load_msa(path):
    msa = []
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
            msa.append((spid, seq))
    return msa


tree_template = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
spids = set([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])

OG_filter = pd.read_table('../../ortho_MSA/OG_filter/out/OG_filter.tsv')
df = pd.read_table('../../ortho_MSA/gap_contrasts/out/total_sums.tsv').merge(OG_filter[['OGid', 'sqidnum']], on='OGid', how='left')  # total_sums.tsv has gnidnum already
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

if not os.path.exists('out/norm1/'):
    os.makedirs('out/norm1/')

head = df.sort_values(by='norm1', ascending=False).head(150)
for i, row in enumerate(head.itertuples()):
    # Load msa and trim terminal insertions
    msa = trim_terminals(load_msa(f'../../ortho_MSA/realign_hmmer/out/{row.OGid}.mfa'))

    # Load decoded states
    posterior = []
    with open(f'../trim_extract/out/{row.OGid}.tsv') as file:
        header = file.readline().split()
        for line in file:
            d = {key: float(value) for key, value in zip(header, line.split())}
            posterior.append(d['2'] + d['3'])
    posterior = np.array(posterior)
    gradient = np.gradient(posterior)

    # Make trim plot
    slices = get_slices(msa, posterior, gradient)
    trims = np.zeros(len(posterior))
    for s in slices:
        trims[s] = 1

    tree = tree_template.shear([seq[0] for seq in msa])
    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    msa = [seq[1].upper() for seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
    plot_msa_lines(msa, [posterior, trims])
    plt.savefig(f'out/norm1/{i}_{row.OGid}.png')
    plt.close()

if not os.path.exists('out/norm2/'):
    os.makedirs('out/norm2/')

head = df.sort_values(by='norm2', ascending=False).head(150)
for i, row in enumerate(head.itertuples()):
    # Load msa and trim terminal insertions
    msa = trim_terminals(load_msa(f'../../ortho_MSA/realign_hmmer/out/{row.OGid}.mfa'))

    # Load decoded states
    posterior = []
    with open(f'../trim_extract/out/{row.OGid}.tsv') as file:
        header = file.readline().split()
        for line in file:
            d = {key: float(value) for key, value in zip(header, line.split())}
            posterior.append(d['2'] + d['3'])
    posterior = np.array(posterior)
    gradient = np.gradient(posterior)

    # Make trim plot
    slices = get_slices(msa, posterior, gradient)
    trims = np.zeros(len(posterior))
    for s in slices:
        trims[s] = 1

    tree = tree_template.shear([seq[0] for seq in msa])
    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    msa = [seq[1].upper() for seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
    plot_msa_lines(msa, [posterior, trims])
    plt.savefig(f'out/norm2/{i}_{row.OGid}.png')
    plt.close()

"""
DEPENDENCIES
../../ortho_MSA/realign_hmmer/realign_hmmer.py
    ../../ortho_MSA/realign_hmmer/out/*
./decode.py
    ./out/*.tsv
"""