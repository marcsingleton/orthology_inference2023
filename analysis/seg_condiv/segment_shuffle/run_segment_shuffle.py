"""Execute the segment_shuffle.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/segment_shuffle.py conserved', shell=True)

"""
DEPENDENCIES
../../../src/segment_shuffle.py
../feature_calc/feature_calc.py
    ../feature_calc/sequences_*.tsv
"""