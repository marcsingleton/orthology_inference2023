"""Execute the segment_shuffle.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/segment_shuffle.py ordered', shell=True)

"""
DEPENDENCIES
../../../src/segment_shuffle.py
../feature_calc/feature_calc.py
    ../feature_calc/sequences_*.tsv
"""