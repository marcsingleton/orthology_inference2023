"""Execute the len_dists.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/len_dists.py ../segment_iupred2a/segment_iupred2a.tsv ordered Ordered Disordered', shell=True)

"""
OUTPUT
Ordered First Ten Counts
1     2720
2     5071
3     6270
4     6681
5     6691
6     6131
7     5547
8     5194
9     4553
10    4142
Name: seq, dtype: int64

Disordered First Ten Counts
1     4026
2     6714
3     8545
4     9345
5     9170
6     8426
7     7453
8     6755
9     5939
10    5140
Name: seq, dtype: int64

NOTES
The ordered and disordered length distributions have heavier tails than their conserved and diverged counterparts
    This is reflected in the greater number of conserved and diverged subsequences; more subsequences results in shorter lengths on average

DEPENDENCIES
../../../src/len_dists.py
../segment_iupred2a/segment_iupred2a.PY
    ../segment_iupred2a/segment_iupred2a.tsv
"""