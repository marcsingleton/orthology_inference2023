"""Execute the aa_pca.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_pca.py ../sample_segs/out/ conserved Conserved Diverged', shell=True)

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
1: [0.13206279 0.12318378 0.09863576 0.09003079 0.07873048]
2: [0.13814809 0.12384513 0.10223027 0.0930317  0.08074934]
4: [0.14685271 0.12365444 0.10176489 0.09546811 0.07962271]
8: [0.16367689 0.12735185 0.11043737 0.10014565 0.0767275 ]
16: [0.17951119 0.12942898 0.12038227 0.09261545 0.07246496]
32: [0.18487645 0.14180928 0.11962467 0.08734077 0.06946414]

NOTES
"Spokes" in PCA are likely short homopolymeric sequences as there are fewer in plots with higher thresholds
At higher cutoffs, the plots are largely the same, though the spokes are more pronounced in the diverged subsequences
    As diverged subsequences are shorter, homopolymeric repeats are likely more prevalent
    Some homopolymeric regions are likely present in all alignments (though the length may differ)

DEPENDENCIES
../../../src/aa_pca.py
../sample_segs/sample_segs.py
    ../segment_iupred2a/out/segments_*.tsv
"""