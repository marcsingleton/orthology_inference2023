"""Perform PCA and visualization of conserved and divergent subsequences."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
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
pcas = []
for i in [2 ** x for x in range(6)]:
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

    # Calculate PCAs and transform subsets
    pca = PCA(n_components=n_components)
    pca.fit(TF_array)
    T_pca = pca.transform(T_array)
    F_pca = pca.transform(F_array)
    pcas.append(pca)  # Store model for later output

    # Plot PCAs
    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.8)
    ax.scatter(T_pca[:, 0], T_pca[:, 1], s=2, alpha=0.1, label=T_name)
    ax.scatter(F_pca[:, 0], F_pca[:, 1], s=2, alpha=0.1, label=F_name)
    ax.set_title(f'PCA of Amino Acid Fractions\nin {T_name} and {F_name} Subsequences')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    leg = fig.legend(bbox_to_anchor=(0.8, 0.5), loc="center left", markerscale=2.5, handletextpad=0)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'aa_pca{i}_combined.png')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.15)
    ax1.scatter(T_pca[:, 0], T_pca[:, 1], s=2, alpha=0.1, label=T_name, color='C0')
    ax2.scatter(F_pca[:, 0], F_pca[:, 1], s=2, alpha=0.1, label=F_name, color='C1')
    ax1.set_title(f'PCA of Amino Acid Fractions\nin {T_name} and {F_name} Subsequences')
    ax2.set_xlabel('PC1')
    ax1.set_ylabel('PC2')
    ax2.set_ylabel('PC2')
    leg = fig.legend(bbox_to_anchor=(0.525, 0), loc="lower center", ncol=2, markerscale=2.5)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'aa_pca{i}_separate.png')

# Print explained variance of components
print(f'Explained variance ratio of first {n_components} components by length cutoff')
for i, pca in enumerate(pcas):
    print(2 ** i, pca.explained_variance_ratio_, sep=': ')
