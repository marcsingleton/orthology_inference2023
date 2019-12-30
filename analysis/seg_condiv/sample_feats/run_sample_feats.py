"""Execute the sample_feats.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/sample_feats.py 3 conserved', shell=True)

"""
DEPENDENCIES
../../src/sample_feats.py
../sample_segs/run_sample_segs.py
    ../sample_segs/segments_*.tsv
"""