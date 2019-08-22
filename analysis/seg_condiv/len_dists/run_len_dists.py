"""Execute the len_dists.py script using the alignment score segmented segmented sequences as input."""

from subprocess import run

run('python ../../../src/len_dists.py ../segment_aliscore/segment_aliscore_ungap.tsv conserved Conserved Diverged', shell=True)

"""
OUTPUT
Conserved First Ten Counts
0    10058
1     3931
2     6401
3     7988
4     9295
5     8360
6     7811
7     7291
8     6526
9     6164
Name: seq, dtype: int64

Diverged First Ten Counts
0    64991
1    31707
2    29264
3    25332
4    21825
5    18023
6    15203
7    12593
8    10864
9     9105
Name: seq, dtype: int64

DEPENDENCIES
../../../src/len_dists.py
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/segment_aliscore_ungap.tsv
"""