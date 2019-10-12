"""Execute the aa_pca.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_pca.py ../segment_aliscore/segment_aliscore_ungap.tsv conserved Conserved Diverged', shell=True)

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
1: [0.1305072  0.12551688 0.10045191 0.09082807 0.07945145]
2: [0.13352252 0.12957343 0.10009022 0.09076264 0.08101013]
4: [0.14725188 0.12965084 0.10409594 0.09934824 0.07996019]
8: [0.16786116 0.12565534 0.11166224 0.09854166 0.07705222]
16: [0.17760895 0.13070081 0.12216715 0.08805334 0.07371074]
32: [0.18459773 0.14147164 0.11737025 0.08688669 0.07027735]

NOTES
"Spokes" in PCA are likely short homopolymeric sequences as there are fewer in plots with higher thresholds
At higher cutoffs, the plots are largely the same, though the spokes are more pronounced in the diverged subsequences
    As diverged subsequences are shorter, homopolymeric repeats are likely more prevalent
    Some homopolymeric regions are likely present in all alignments (though the length may differ)

DEPENDENCIES
../../../src/aa_pca.py
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/segment_aliscore_ungap.tsv
"""