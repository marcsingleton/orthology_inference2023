"""Execute the aa_pca.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_pca.py ../segment_iupred2a/segment_iupred2a.tsv ordered Ordered Disordered', shell=True)

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
1: [0.11084521 0.09524703 0.07807074 0.07550875 0.06941813]
2: [0.1184349  0.10076891 0.08082267 0.07834924 0.07228805]
4: [0.12314819 0.10819302 0.08149624 0.08029265 0.07385251]
8: [0.13346104 0.105507   0.09238172 0.08444683 0.07533771]
16: [0.14505671 0.10597107 0.09950995 0.0843975  0.07703217]
32: [0.16456761 0.12873541 0.09513586 0.08929247 0.07620161]

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