"""Execute the feature_calc.py script using the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_calc_scripts/feature_calc.py ../segment_aliscore/out/segment_aliscore.tsv conserved', shell=True)

"""
NOTES
There are 8560 diverged sequences with len >= 32.

DEPENDENCIES
../../../src/feature_calc/feature_calc.py
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/out/segment_aliscore.tsv
"""