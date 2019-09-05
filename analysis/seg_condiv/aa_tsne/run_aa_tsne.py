"""Execute the aa_tsne.py script using the alignment score segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_tsne.py ../segment_aliscore/segment_aliscore_ungap.tsv conserved Conserved Diverged', shell=True)

"""
DEPENDENCIES
../../../src/aa_tsne.py
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/segment_aliscore_ungap.tsv
"""