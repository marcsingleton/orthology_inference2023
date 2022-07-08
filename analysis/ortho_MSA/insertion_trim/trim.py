"""Trim state 1 segments to yield trimmed alignments."""

import numpy as np
from src.ortho_MSA.trim import get_slices
from src.utils import read_fasta

posterior_high = 0.75
posterior_low = 0.5
gradient_high = 0.02
gradient_low = 0.001

# Load OGids
OGids = []
with open('../realign_hmmer/out/errors.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGid, error_flag = fields['OGid'], fields['error_flag']
        if error_flag == 'False':
            OGids.append(OGid)

for OGid in OGids:
    msa = read_fasta(f'../realign_hmmer/out/mafft/{OGid}.afa')

    # Load decoded states and calculate derivative
    posterior = []
    with open(f'out/{OGid}.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: float(value) for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            posterior.append(fields['2'] + fields['3'])
    posterior = np.array(posterior)
    gradient = np.gradient(posterior)

    # Find trimmed regions
    slices = get_slices(msa, posterior, gradient, posterior_high, posterior_low, gradient_high, gradient_low)

    # Invert slices
    # Original slices may overlap due to extension,
    # but this will produce empty "inverted" slices and will not affect the result.
    inverts, i = [], 0
    for s in slices:
        inverts.append(slice(i, s.start))
        i = s.stop
    inverts.append(slice(i, len(msa[0][1])))

    # Write trimmed MSA
    with open(f'out/{OGid}.afa', 'w') as file:
        for header, seq1 in msa:
            seq2 = ''.join([seq1[s] for s in inverts])
            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)])
            file.write(f'{header}\n{seqstring}\n')

"""
DEPENDENCIES
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/errors.tsv
    ../realign_hmmer/out/mafft/*.afa
./decode.py
    ./out/*.tsv
"""