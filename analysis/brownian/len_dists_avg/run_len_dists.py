"""Execute the len_dists.py script using the IUPRED2a average segmented sequences as input."""

from subprocess import run

run('python ../../../src/len_dists.py ../segment_avg/segment_avg.tsv ordered Ordered Disordered', shell=True)

"""
OUTPUT
Ordered First Ten Counts
0    20760
1     4118
2     5636
3     6614
4     6671
5     6452
6     6354
7     5346
8     4991
9     4420
Name: seq, dtype: int64

Disordered First Ten Counts
0    20944
1     5414
2     8450
3     8917
4     9558
5     8761
6     8357
7     7315
8     6521
9     5705
Name: seq, dtype: int64

DEPENDENCIES
../../../src/len_dists.py
../segment_avg/segment_avg.py
    ../segment_avg/segment_avg.tsv
"""