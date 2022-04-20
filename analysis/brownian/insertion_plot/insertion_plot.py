"""Plot decoded insertions in alignments with large gap contrasts."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio
from src.draw import plot_msa_data
from src.brownian2.trim import trim_terminals, get_slices
from src.utils import read_fasta

tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}
spids = {tip.name for tip in tree.tips() if tip.name != 'sleb'}

df = pd.read_table('../../ortho_MSA/gap_contrasts/out/total_sums.tsv')
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

for label in ['norm1', 'norm2']:
    if not os.path.exists(f'out/{label}/'):
        os.makedirs(f'out/{label}/')

    head = df.sort_values(by=label, ascending=False).head(150)
    for i, row in enumerate(head.itertuples()):
        # Load msa and trim terminal insertions
        msa = trim_terminals(read_fasta(f'../../ortho_MSA/realign_hmmer/out/{row.OGid}.afa'))
        msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]

        # Load decoded states
        posterior = []
        with open(f'../insertion_trim/out/{row.OGid}.tsv') as file:
            header = file.readline().rstrip('\n').split('\t')
            for line in file:
                fields = {key: float(value) for key, value in zip(header, line.rstrip('\n').split('\t'))}
                posterior.append(fields['2'] + fields['3'])
        posterior = np.array(posterior)
        gradient = np.gradient(posterior)

        # Make trim plot
        slices = get_slices(msa, posterior, gradient)
        trims = np.zeros(len(posterior))
        for s in slices:
            trims[s] = 1

        msa = [seq.upper() for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only
        plot_msa_data(msa, [posterior, trims], figsize=(15, 6))
        plt.savefig(f'out/{label}/{i:03}_{row.OGid}.png', bbox_inches='tight')
        plt.close()

"""
DEPENDENCIES
../../ortho_MSA/gap_contrasts/gap_contrasts.py
    ../gap_contrasts/out/total_sums.tsv
../../ortho_MSA/realign_hmmer/realign_hmmer.py
    ../../ortho_MSA/realign_hmmer/out/*
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../insertion_trim/decode.py
    ./insertion_trim/out/*.tsv
"""