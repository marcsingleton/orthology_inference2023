"""Execute the aa_tsne.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_tsne.py ../sample_segs/out/ ordered Ordered Disordered', shell=True)

"""
DEPENDENCIES
../../../src/aa_tsne.py
../sample_segs/sample_segs.py
    ../segment_iupred2a/out/segments_*.tsv
"""