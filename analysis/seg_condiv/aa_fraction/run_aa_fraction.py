"""Execute the aa_fraction.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_fraction.py ../segment_aliscore/segment_aliscore.tsv conserved Conserved Diverged', shell=True)

"""
OUTPUT
Number of conserved subsequences: 322940
Number of diverged subsequences: 305170
Number of conserved amino acids: 25887032
Number of diverged amino acids: 2197778

DEPENDENCIES
../../../src/aa_fraction.py
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/segment_aliscore_ungap.tsv
"""