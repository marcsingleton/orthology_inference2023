"""Execute the aa_fraction.py script using the IUPRED2a average segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_fraction.py ../segment_avg/segment_avg.tsv ordered Ordered Disordered', shell=True)

"""
OUTPUT
Number of ordered subsequences: 242440
Number of disordered subsequences: 221470
Number of ordered amino acids: 20275996
Number of disordered amino acids: 7808814

DEPENDENCIES
../../../src/aa_fraction.py
../segment_iupred2a/segment_avg.py
    ../segment_iupred2a/segment_avg.tsv
"""