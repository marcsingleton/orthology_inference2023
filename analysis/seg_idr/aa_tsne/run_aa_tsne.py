"""Execute the aa_tsne.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/aa_tsne.py ../segment_iupred2a/segment_iupred2a.tsv ordered Ordered Disordered', shell=True)

"""
DEPENDENCIES
../../../src/aa_tsne.py
../segment_iupred2a/segment_iupred2a.py
    ../segment_iupred2a/segment_iupred2a.tsv
"""