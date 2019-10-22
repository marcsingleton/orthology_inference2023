"""Perform t-SNE and visualization of two types of subsequences."""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import re
from sklearn.manifold import TSNE
from sys import argv


def fracs(seq):
    return {char: seq.count(char) / len(seq) for char in alphabet}


# Input variables
segment_dir = argv[1]  # Segment directory must end in /
type_name = argv[2]  # Name of column denoting segment type
T_name = argv[3]  # Name of True type in sentence case
F_name = argv[4]  # Name of False type in sentence case

# Constants
alphabet = 'DEHKRNQSTAILMVFWYCGP'
n_components = 5

paths = filter(lambda x: re.match('segments_[0-9]+\.tsv', x), os.listdir(segment_dir))
for path in paths:
    kl_divs = []

    # Read data and split segments
    segs = pd.read_csv(segment_dir + path, sep='\t', keep_default_na=False)
    T_seqs = segs.loc[segs[type_name], 'seq']
    F_seqs = segs.loc[~segs[type_name], 'seq']  # ~ is bitwise NOT operator; it interacts properly with numpy objects but not Python booleans

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = path[j0 + 1:j1]

    # Calculate fractions and convert to array
    T_features = T_seqs.apply(fracs)
    F_features = F_seqs.apply(fracs)
    T_array = np.array([[*entry.values()] for entry in T_features])
    F_array = np.array([[*entry.values()] for entry in F_features])
    array = np.concatenate((T_array, F_array), axis=0)

    # Make output directories for feature sets
    cur_dir = f'{i}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    for j in [5 * 2 ** x for x in range(8)]:
        # Calculate PCAs and transform subsets
        tsne = TSNE(n_components=2, perplexity=j)
        transform = tsne.fit_transform(array)
        T_tsne = transform[segs[type_name], :]
        F_tsne = transform[~segs[type_name], :]
        kl_divs.append((j, tsne.kl_divergence_))  # Store model for later output

        # Plot t-SNEs
        # One panel
        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.175)
        ax.scatter(T_tsne[:, 0], T_tsne[:, 1], s=2, alpha=0.1, label=T_name)
        ax.scatter(F_tsne[:, 0], F_tsne[:, 1], s=2, alpha=0.1, label=F_name)
        ax.set_title(f't-SNE of Amino Acid Fractions\nin {T_name} and {F_name} Subsequences')
        leg = fig.legend(bbox_to_anchor=(0.5, 0), loc='lower center', ncol=2, markerscale=2.5)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        plt.savefig(cur_dir + f'aa_tsne{j}_combined.png')
        plt.close()

        # Two panels
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.15)
        ax1.scatter(T_tsne[:, 0], T_tsne[:, 1], s=2, alpha=0.1, label=T_name, color='C0')
        ax2.scatter(F_tsne[:, 0], F_tsne[:, 1], s=2, alpha=0.1, label=F_name, color='C1')
        ax1.set_title(f't-SNE of Amino Acid Fractions\nin {T_name} and {F_name} Subsequences')
        leg = fig.legend(bbox_to_anchor=(0.525, 0), loc='lower center', ncol=2, markerscale=2.5)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        plt.savefig(cur_dir + f'aa_tsne{j}_separate.png')
        plt.close()

        # Write model
        with open(cur_dir + f'model{j}.pickle', 'wb') as file:
            pickle.dump(tsne, file)

    # Write model summaries
    with open(cur_dir + 'model_summary.txt', 'w') as file:
        file.write('#perplexity\tKL_divergence\n')
        for perp, kl_div in kl_divs:
            file.write(f'{perp}\t{kl_div}\n')
