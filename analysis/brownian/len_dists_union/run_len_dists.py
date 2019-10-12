"""Execute the len_dists.py script using the IUPRED2a union segmented sequences as input."""

from subprocess import run

run('python ../../../src/len_dists.py ../segment_union/segment_union.tsv ordered Ordered Disordered', shell=True)

"""
OUTPUT
Ordered First Ten Counts
0    29207
1    17994
2    16583
3    15158
4    13698
5    12181
6     9806
7     8861
8     7638
9     7088
Name: seq, dtype: int64

Disordered First Ten Counts
0    87544
1    19692
2    16428
3    13981
4    13389
5    11237
6     9968
7     8237
8     7725
9     6808

NOTES
Compared to segmentation by average, segmentation by union creates more short sequences, particularly for disordered regions.
    This is unsurprising for disordered regions since they can be disrupted by a single ordered subsequence.
    The increase in short ordered regions suggests the presence of disordered insertions in a single sequence.

DEPENDENCIES
../../../src/len_dists.py
../segment_union/segment_union.py
    ../segment_union/segment_union.tsv
"""