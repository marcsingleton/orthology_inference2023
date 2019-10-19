"""Execute the sample_segs.py script using the IUPRED2a segmented sequences as input."""

from subprocess import run

run('python ../../../src/sample_segs.py ../segment_iupred2a/segment_iupred2a.tsv', shell=True)

"""
DEPENDENCIES
../../../src/sample_segs.py
../segment_iupred2a/segment_iupred2a.py
    ../segment_iupred2a/segment_iupred2a.tsv
"""