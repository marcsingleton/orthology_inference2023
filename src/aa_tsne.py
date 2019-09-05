"""Perform t-SNE and visualization of conserved and divergent subsequences."""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
from sklearn.manifold import TSNE
from sys import argv


def fracs(seq):
    return {char: seq.count(char) / len(seq) for char in alphabet}


# Input variables
path = argv[1]  # Path to segmented sequences .tsv
key = argv[2]  # Key of column denoting subsequences class
T_name = argv[3]  # Name of True class in sentence case
F_name = argv[4]  # Name of False class in sentence case

# Constants
alphabet = 'DEHKRNQSTAILMVFWYCGP'
n_components = 5

df = pd.read_csv(path, sep='\t', keep_default_na=False)
for i in [2 ** x for x in range(6)]:
    kl_divs = []

    # Read data and split subsequences
    df = df[df['seq'].map(lambda x: len(x) >= i)]  # Select entries where sequence is above length threshold
    df_T = df[df[key] == True]
    df_F = df[df[key] == False]

    # Sample full dataframes, calculate fractions, and convert to array
    T_sample = df_T.sample(8500)['seq']
    F_sample = df_F.sample(8500)['seq']
    T_features = T_sample.apply(fracs)
    F_features = F_sample.apply(fracs)
    T_array = np.array([[*entry.values()] for entry in T_features])
    F_array = np.array([[*entry.values()] for entry in F_features])
    TF_array = np.concatenate((T_array, F_array), axis=0)

    # Make output directories for feature sets
    cur_dir = f'{i}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    for j in [5 * 2 ** x for x in range(8)]:
        # Calculate PCAs and transform subsets
        tsne = TSNE(n_components=2, perplexity=j)
        trans = tsne.fit_transform(TF_array)
        T_tsne = trans[:8500, :]
        F_tsne = trans[8500:, :]
        kl_divs.append((j, tsne.kl_divergence_))  # Store model for later output

        # Plot t-SNEs
        fig, ax = plt.subplots()
        fig.subplots_adjust(right=0.8)
        ax.scatter(T_tsne[:, 0], T_tsne[:, 1], s=2, alpha=0.1, label=T_name)
        ax.scatter(F_tsne[:, 0], F_tsne[:, 1], s=2, alpha=0.1, label=F_name)
        ax.set_title(f't-SNE of Amino Acid Fractions\nin {T_name} and {F_name} Subsequences')
        leg = fig.legend(bbox_to_anchor=(0.8, 0.5), loc="center left", markerscale=2.5, handletextpad=0)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        plt.savefig(cur_dir + f'aa_tsne{j}_combined.png')
        plt.close()

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.15)
        ax1.scatter(T_tsne[:, 0], T_tsne[:, 1], s=2, alpha=0.1, label=T_name, color='C0')
        ax2.scatter(F_tsne[:, 0], F_tsne[:, 1], s=2, alpha=0.1, label=F_name, color='C1')
        ax1.set_title(f't-SNE of Amino Acid Fractions\nin {T_name} and {F_name} Subsequences')
        leg = fig.legend(bbox_to_anchor=(0.525, 0), loc="lower center", ncol=2, markerscale=2.5)
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
