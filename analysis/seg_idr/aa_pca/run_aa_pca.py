"""Execute the aa_pca.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_pca.py ../segment_iupred2a/segment_iupred2a.tsv ordered Ordered Disordered', shell=True)

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
1: [0.11202364 0.10079852 0.0787125  0.07605646 0.06953976]
2: [0.1138936  0.10402136 0.08039153 0.0769274  0.07168185]
4: [0.12214567 0.1032966  0.08654573 0.07993077 0.07456813]
8: [0.13283162 0.10291754 0.09406589 0.08461914 0.0750315 ]
16: [0.15037903 0.10966678 0.09557253 0.08403114 0.07549575]
32: [0.16305408 0.12818101 0.09469934 0.09164521 0.07413461]

NOTES
Spokes are much less prevalent in the IDR PCA than in the conserved/diverged PCA
    This is likely an effect of the average longer subsequence lengths in the IDR set; longer lengths results in fewer di or tri peptide sequences
At higher length cutoffs, the disordered sequences appear to slightly separate from the ordered sequences
    Perhaps at higher cutoffs, their intrinsic amino acid preferences begin to overcome the noise inherent in small sequences

DEPENDENCIES
../../../src/aa_pca.py
../segment_iupred2a/segment_iupred2a.py
    ../segment_iupred2a/segment_iupred2a.tsv
"""