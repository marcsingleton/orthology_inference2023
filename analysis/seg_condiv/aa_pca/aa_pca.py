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
type_name = 'conserved'  # Name of column denoting segment type
T_name = 'Conserved'  # Name of True type in sentence case
F_name = 'Diverged'  # Name of False type in sentence case

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
1: [0.1284349  0.12286042 0.09960856 0.09187439 0.07837038]
2: [0.13269481 0.12585695 0.10107725 0.09484571 0.07881363]
4: [0.15472376 0.12399268 0.1058095  0.09464213 0.07670636]
8: [0.16578987 0.13025587 0.10965513 0.09982943 0.07513263]
16: [0.19262305 0.1274957  0.11882532 0.08707529 0.07050748]
32: [0.1882316  0.13956229 0.11619777 0.0860077  0.06807423]

NOTES
"Spokes" in PCA are likely short homopolymeric sequences as there are fewer in plots with higher thresholds
At higher cutoffs, the plots are largely the same, though the spokes are more pronounced in the diverged subsequences
    As diverged subsequences are shorter, homopolymeric repeats are likely more prevalent
    Some homopolymeric regions are likely present in all alignments (though the length may differ)

DEPENDENCIES
../sample_segs/sample_segs.py
    ../segment_iupred2a/out/segments_*.tsv
"""