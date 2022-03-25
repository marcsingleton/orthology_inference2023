"""Extract state 1 segments to yield trimmed alignments."""

import os

import numpy as np
from src.brownian2.trim import trim_terminals, get_slices
from src.utils import read_fasta


OGids = [path.split('.')[0] for path in os.listdir('../../ortho_MSA/realign_hmmer/out/') if path.endswith('.mfa')]
for OGid in OGids:
    # Load msa and trim terminal insertions
    msa = trim_terminals(read_fasta(f'../../ortho_MSA/realign_hmmer/out/{OGid}.afa'))

    # Load decoded states and calculate derivative
    posterior = []
    with open(f'out/{OGid}.tsv') as file:
        header = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: float(value) for key, value in zip(header, line.rstrip('\n').split('\t'))}
            posterior.append(fields['2'] + fields['3'])
    posterior = np.array(posterior)
    gradient = np.gradient(posterior)

    # Find trimmed regions
    slices = get_slices(msa, posterior, gradient)

    # Invert slices
    # Original slices may overlap due to extension,
    # but this will produce empty "inverted" slices and will not affect the result.
    inverts, i = [], 0
    for s in slices:
        inverts.append(slice(i, s.start))
        i = s.stop
    inverts.append(slice(i, len(msa[0][1])))

    # Write trimmed MSA
    with open(f'out/{OGid}.mfa', 'w') as file:
        for header, seq1 in msa:
            seq2 = ''.join([seq1[s] for s in inverts])
            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)]) + '\n'
            file.write(header + '\n' + seqstring)

"""
DEPENDENCIES
../../ortho_MSA/realign_hmmer/realign_hmmer.py
    ../../ortho_MSA/realign_hmmer/out/*.afa
./decode.py
    ./out/*.tsv
"""