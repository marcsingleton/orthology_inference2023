"""Perform PCA and visualization of conserved and divergent subsequences."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def fracs(seq):
    return {char: seq.count(char) / len(seq) for char in alphabet}


path = '../segment_aliscore/segment_aliscore_ungap.tsv'
alphabet = 'DEHKRNQSTAILMVFWYCGP'
n_components = 5


df = pd.read_csv(path, sep='\t', keep_default_na=False)
pcas = []
for i in range(10):
    # Split into separate conserved and divergent dataframes
    df = df[df['seq'].map(lambda x: len(x) > i)]  # Select entries where sequence is above length threshold
    df_con = df[df['conserved'] == True]
    df_div = df[df['conserved'] == False]

    # Sample full dataframes, calculate fractions, and convert to array
    con_sample = df_con.sample(10000)['seq']
    div_sample = df_div.sample(10000)['seq']
    con_features = con_sample.apply(fracs)
    div_features = div_sample.apply(fracs)
    con_array = np.array([[*entry.values()] for entry in con_features])
    div_array = np.array([[*entry.values()] for entry in div_features])
    condiv_array = np.concatenate((con_array, div_array), axis=0)

    # Calculate PCAs and transform subsets
    pca = PCA(n_components=n_components)
    pca.fit(condiv_array)
    con_pca = pca.transform(con_array)
    div_pca = pca.transform(div_array)
    pcas.append(pca)  # Store model for later output

    # Plot PCAs
    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.8)
    ax.scatter(con_pca[:, 0], con_pca[:, 1], s=2, alpha=0.1, label='Conserved')
    ax.scatter(div_pca[:, 0], div_pca[:, 1], s=2, alpha=0.1, label='Diverged')
    ax.set_title('PCA of Amino Acid Fractions\nin Conserved and Diverged Subsequences')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    leg = fig.legend(bbox_to_anchor=(0.8, 0.5), loc="center left", markerscale=2.5, handletextpad=0)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.savefig(f'aa_pca{i}_combined.png')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.15)
    ax1.scatter(con_pca[:, 0], con_pca[:, 1], s=2, alpha=0.1, label='Conserved', color='C0')
    ax2.scatter(div_pca[:, 0], div_pca[:, 1], s=2, alpha=0.1, label='Diverged', color='C1')
    ax1.set_title('PCA of Amino Acid Fractions\nin Conserved and Diverged Subsequences')
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
    print(i, pca.explained_variance_ratio_, sep=': ')

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
0: [0.12815369 0.12236204 0.10081884 0.0911327  0.07762422]
1: [0.13668873 0.12634128 0.1017678  0.09461824 0.08159819]
2: [0.14253811 0.12855066 0.10286692 0.09594496 0.07928601]
3: [0.14895979 0.12751305 0.1036118  0.09631596 0.07874296]
4: [0.14789266 0.1281257  0.10498759 0.09897445 0.07940961]
5: [0.1590381  0.12537549 0.10867264 0.09741331 0.07748833]
6: [0.16128051 0.12915957 0.1105486  0.09912541 0.07672047]
7: [0.16375079 0.1275373  0.11232373 0.09893851 0.07807543]
8: [0.16736497 0.12776024 0.11385524 0.09826462 0.07513673]
9: [0.17393421 0.12895266 0.11593678 0.09184615 0.07625517]

NOTES
"Spokes" in PCA are likely short homopolymeric sequences as there are fewer in plots with higher thresholds
At higher cutoffs, the plots are largely the same, though the spokes are more pronounced in the diverged subsequences
    As diverged subsequences are shorter, homopolymeric repeats are likely more prevalent
    Some homopolymeric regions are likely present in all alignments (though the length may differ)

DEPENDENCIES
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/segment_aliscore_ungap.tsv
"""