"""Execute the segment_shuffle.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/segment_shuffle.py ../sample_segs/out/ conserved', shell=True)

"""
DEPENDENCIES
../../../src/segment_shuffle.py
../sample_segs/sample_seg.py
    ../sample_segs/out/segments_*.tsv
"""