"""Execute the aa_pca.py script using the IUPRED2a average segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_pca.py ../segment_avg/segment_avg.tsv ordered Ordered Disordered', shell=True)

"""
OUTPUT
Explained variance ratio of first 5 components by length cutoff
1: [0.09665858 0.08828238 0.07974026 0.0715485  0.06560808]
2: [0.10634983 0.09372372 0.08320178 0.0714325  0.06733368]
4: [0.10743563 0.09262671 0.09102385 0.07386548 0.07098653]
8: [0.11818862 0.1101219  0.09560506 0.07554226 0.07034325]
16: [0.12496676 0.12067177 0.0941604  0.08654011 0.07250974]
32: [0.1475941  0.13407494 0.10013184 0.08233643 0.07747892]

DEPENDENCIES
../../../src/aa_pca.py
../segment_iupred2a/segment_avg.py
    ../segment_iupred2a/segment_avg.tsv
"""