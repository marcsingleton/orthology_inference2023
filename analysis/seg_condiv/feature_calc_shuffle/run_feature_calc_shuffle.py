"""Execute the feature_calc_shuffle.py script using the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_calc_scripts/feature_calc_shuffle.py ../segment_shuffle/out/ conserved', shell=True)

"""
DEPENDENCIES
../../../src/feature_calc/feature_calc_shuffle.py
../segment_shuffle/segment_shuffle.py
    ../segment_shuffle/out/shuffseq_*.tsv
"""