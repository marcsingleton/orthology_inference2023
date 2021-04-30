"""Execute the dists.py script using the features of the unshuffled sequences as input."""

from subprocess import run

run('python ../../../src/seg_scripts/plot/dists.py ../sample_feats/out/ 3 ordered Ordered Disordered', shell=True)

"""
DEPENDENCIES
../../../src/seg_scripts/plot/dists.py
../sample_feats/run_sample_feats.py
    ../sample_feats/out/features_*.tsv
"""