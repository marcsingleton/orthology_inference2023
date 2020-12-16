"""Perform PCA and visualization of two types of subsequences."""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
from sklearn.decomposition import PCA


def fracs(seq):
    return {char: seq.count(char) / len(seq) for char in alphabet}


# Input variables
segment_dir = '../sample_segs/out/'  # Segment directory must end in /
type_name = 'ordered'  # Name of column denoting segment type
T_name = 'Ordered'  # Name of True type in sentence case
F_name = 'Disordered'  # Name of False type in sentence case

# Constants
alphabet = 'DEHKRNQSTAILMVFWYCGP'
n_components = 5

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

pcas = {}
paths = filter(lambda x: re.match('segments_[0-9]+\.tsv', x), os.listdir(segment_dir))
for path in paths:
    # Read data and split segments
    segs = pd.read_csv(segment_dir + path, sep='\t', keep_default_na=False)
    T_seqs = segs.loc[segs[type_name], 'seq']
    F_seqs = segs.loc[~segs[type_name], 'seq']  # ~ is bitwise NOT operator; it interacts properly with numpy objects but not Python booleans

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = int(path[j0 + 1:j1])

    # Calculate fractions and convert to array
    T_features = T_seqs.apply(fracs)
    F_features = F_seqs.apply(fracs)
    T_array = np.array([[*entry.values()] for entry in T_features])
    F_array = np.array([[*entry.values()] for entry in F_features])
    array = np.concatenate((T_array, F_array), axis=0)

    # Calculate PCAs and transform subsets
    pca = PCA(n_components=n_components)
    pca.fit(array)
    T_pca = pca.transform(T_array)
    F_pca = pca.transform(F_array)
    pcas[i] = pca  # Store model for later output

    # Plot PCAs
    # One panel
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.175)
    ax.scatter(T_pca[:, 0], T_pca[:, 1], s=2, alpha=0.1, label=T_name)
    ax.scatter(F_pca[:, 0], F_pca[:, 1], s=2, alpha=0.1, label=F_name)
    ax.set_title(f'PCA of Amino Acid Fractions\nin {T_name} and {F_name} Subsequences')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    leg = fig.legend(bbox_to_anchor=(0.5, 0), loc='lower center', ncol=2, markerscale=2.5)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'out/aa_pca{i}_combined.png')

    # Two panels
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.15)
    ax1.scatter(T_pca[:, 0], T_pca[:, 1], s=2, alpha=0.1, label=T_name, color='C0')
    ax2.scatter(F_pca[:, 0], F_pca[:, 1], s=2, alpha=0.1, label=F_name, color='C1')
    ax1.set_title(f'PCA of Amino Acid Fractions\nin {T_name} and {F_name} Subsequences')
    ax2.set_xlabel('PC1')
    ax1.set_ylabel('PC2')
    ax2.set_ylabel('PC2')
    leg = fig.legend(bbox_to_anchor=(0.525, 0), loc='lower center', ncol=2, markerscale=2.5)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'out/aa_pca{i}_separate.png')

# Print explained variance of components
print(f'Explained variance ratio of first {n_components} components by length cutoff')
for i, pca in sorted(pcas.items()):
    print(i, pca.explained_variance_ratio_, sep=': ')

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
1: [0.11051914 0.0970163  0.07686716 0.07622957 0.07247313]
2: [0.11356423 0.10164706 0.08359128 0.0787637  0.0717176 ]
4: [0.12260539 0.10260207 0.08730236 0.08021261 0.07131195]
8: [0.13289551 0.10234792 0.09325454 0.0820337  0.07541747]
16: [0.14425587 0.11004906 0.09992752 0.08386489 0.07628767]
32: [0.1644331  0.12567238 0.09494454 0.08828389 0.07491024]

NOTES
Spokes are much less prevalent in the IDR PCA than in the conserved/diverged PCA
    This is likely an effect of the average longer subsequence lengths in the IDR set; longer lengths results in fewer di or tri peptide sequences
At higher length cutoffs, the disordered sequences appear to slightly separate from the ordered sequences
    Perhaps at higher cutoffs, their intrinsic amino acid preferences begin to overcome the noise inherent in small sequences

DEPENDENCIES
../sample_segs/sample_segs.py
    ../segment_iupred2a/out/segments_*.tsv
"""