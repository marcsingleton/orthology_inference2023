"""Execute the sample_segs.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/sample_segs.py ../segment_aliscore/segment_aliscore.tsv', shell=True)

"""
DEPENDENCIES
../../../src/sample_segs.py
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/segment_aliscore.tsv
"""