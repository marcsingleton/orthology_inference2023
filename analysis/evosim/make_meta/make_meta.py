"""Make meta alignments of ordered and disordered regions."""

import os
from collections import namedtuple
from random import randrange, seed


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


def is_redundant(col, cutoff):
    gapnum = 0
    for _, sym in col:
        if sym in ['-', '.', 'X']:
            gapnum += 1
    return gapnum <= (1 - cutoff) * len(col)


Column = namedtuple('Column', ['spid', 'sym'])
seed(930715)  # Set seed to make results consistent

# Load regions
OGid2regions = {}
with open('../../brownian2/aucpred_regions/out/regions.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder = line.split()
        try:
            OGid2regions[OGid].append((int(start), int(stop), disorder))
        except KeyError:
            OGid2regions[OGid] = [(int(start), int(stop), disorder)]

# Extract column pools
colpools = [('100red_D', lambda col: is_redundant(col, 1), []),
            ('100red_O', lambda col: is_redundant(col, 1), []),
            ('50red_D', lambda col: is_redundant(col, 0.5), []),
            ('50red_O', lambda col: is_redundant(col, 0.5), []),
            ('0red_D', lambda col: is_redundant(col, 0), []),
            ('0red_O', lambda col: is_redundant(col, 0), [])]
for OGid, regions in OGid2regions.items():
    msa = load_msa(f'../../brownian2/insertion_trim/out/{OGid}.mfa')
    if len(msa) < 31:  # Only use alignments with all species
        continue

    for start, stop, disorder in regions:
        for i in range(start, stop):
            col = [Column(header[-4:], seq[i]) for header, seq in msa]
            tag = 'D' if disorder == 'True' else 'O'
            for label, condition, colpool in colpools:
                if label[-1] == tag and condition(col):
                    colpool.append(col)

# Make meta alignments
if not os.path.exists('out/'):
    os.mkdir('out/')

for label, _, colpool in colpools:
    print(f'{label}:', len(colpool))
    sample = [colpool[randrange(len(colpool))] for _ in range(int(1E5))]
    seqs = {}
    for col in sample:
        for spid, sym in col:
            try:
                seqs[spid].append(sym)
            except KeyError:
                seqs[spid] = [sym]

    with open(f'out/{label}.fasta', 'w') as file:
        for spid, seq in sorted(seqs.items()):
            seqstring = '\n'.join([''.join(seq[i:i+80]) for i in range(0, len(seq), 80)]) + '\n'
            file.write(f'>{spid} {label}\n' + seqstring)

"""
OUTPUT
100red_D: 655952
100red_O: 3004444
50red_D: 1080267
50red_O: 3294780
0red_D: 1559840
0red_O: 3420469

DEPENDENCIES
../aucpred_regions/aucpred_regions.py
    ../aucpred_regions/out/regions.tsv
../insertion_trim/insertion_trim.py
    ../insertion_trim/out/*.mfa
"""