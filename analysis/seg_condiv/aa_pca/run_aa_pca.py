"""Execute the aa_pca.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_pca.py ../segment_aliscore/segment_aliscore_ungap.tsv conserved Conserved Diverged', shell=True)

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
1: [0.1301792  0.12335692 0.10036955 0.09075553 0.07869939]
2: [0.13423173 0.12724973 0.10090771 0.09635037 0.07897103]
4: [0.14639471 0.12509654 0.10444533 0.09872864 0.07945201]
8: [0.16445647 0.12821788 0.11156149 0.0977362  0.07756936]
16: [0.18209271 0.13235243 0.11977519 0.08843416 0.0726299 ]
32: [0.18437282 0.14026811 0.1192075  0.08803145 0.07037384]

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