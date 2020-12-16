"""Execute the aa_pca.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_pca.py ../sample_segs/out/ ordered Ordered Disordered', shell=True)

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
1: [0.10835237 0.0985546  0.07720514 0.07521693 0.07371808]
2: [0.11563894 0.09935632 0.07967057 0.07658405 0.07078912]
4: [0.12177269 0.10160838 0.08221755 0.0814585  0.07623621]
8: [0.13419047 0.1029753  0.09404481 0.08139625 0.07575718]
16: [0.14722338 0.10811996 0.09972444 0.08523588 0.07748651]
32: [0.16252439 0.12638589 0.09661872 0.08927426 0.07550576]

NOTES
Spokes are much less prevalent in the IDR PCA than in the conserved/diverged PCA
    This is likely an effect of the average longer subsequence lengths in the IDR set; longer lengths results in fewer di or tri peptide sequences
At higher length cutoffs, the disordered sequences appear to slightly separate from the ordered sequences
    Perhaps at higher cutoffs, their intrinsic amino acid preferences begin to overcome the noise inherent in small sequences

DEPENDENCIES
../../../src/aa_pca.py
../sample_segs/sample_segs.py
    ../segment_iupred2a/out/segments_*.tsv
"""