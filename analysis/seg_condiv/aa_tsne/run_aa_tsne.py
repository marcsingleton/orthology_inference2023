"""Execute the aa_tsne.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_tsne.py ../sample_segs/out/ conserved Conserved Diverged', shell=True)

"""
DEPENDENCIES
../../../src/aa_tsne.py
../sample_segs/sample_segs.py
    ../segment_iupred2a/out/segments_*.tsv
"""