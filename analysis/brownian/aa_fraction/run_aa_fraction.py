"""Execute the aa_fraction.py script using the segments corresponding to the pics as input."""

from subprocess import run

run('python ../../../src/aa_fraction.py ../filter_blocks/seg_filter.tsv ordered Ordered Disordered', shell=True)

"""
OUTPUT
Number of ordered subsequences: 196490
Number of disordered subsequences: 173670
Number of ordered amino acids: 19593513
Number of disordered amino acids: 7098410

DEPENDENCIES
../../../src/aa_fraction.py
../filter_blocks/filter_blocks.py
    ../filter_blocks/seg_filter.tsv
"""