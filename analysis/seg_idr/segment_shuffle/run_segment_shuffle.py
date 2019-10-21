"""Execute the segment_shuffle.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/segment_shuffle.py ../sample_segs/ ordered', shell=True)

"""
DEPENDENCIES
../../../src/segment_shuffle.py
../sample_segs/sample_seg.py
    ../sample_segs/segments_*.tsv
"""