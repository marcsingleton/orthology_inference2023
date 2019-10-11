"""Execute the feature_calc.py script using the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/feature_calc_scripts/feature_calc.py ../segment_iupred2a/segment_iupred2a.tsv ordered', shell=True)

"""
DEPENDENCIES
../../../src/feature_calc/feature_calc.py
../segment_iupred2a/segment_iupred2a.py
    ../segment_iupred2a/segment_iupred2a.tsv
"""