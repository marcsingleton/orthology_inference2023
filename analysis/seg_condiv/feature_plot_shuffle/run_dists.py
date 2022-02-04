"""Execute the dists.py script using the features of the shuffled sequences as input."""

from subprocess import run

run('python ../../../src/seg_scripts/plot/dists.py ../feature_calc_shuffle/out/ 1 conserved Conserved Diverged', shell=True)

"""
DEPENDENCIES
../feature_calc_shuffle/run_feature_calc_shuffle.py
    ../feature_calc_shuffle/out/features_*.tsv
"""