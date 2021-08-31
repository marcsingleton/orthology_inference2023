"""Extract state 1 segments to yield trimmed alignments."""

import os

import numpy as np
from src.brownian2.trim import trim_terminals, get_slices


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


OGids = [path.split('.')[0] for path in os.listdir('../../ortho_MSA/realign_hmmer/out/') if path.endswith('.mfa')]
for OGid in OGids:
    # Load msa and trim terminal insertions
    msa = trim_terminals(load_msa(f'../../ortho_MSA/realign_hmmer/out/{OGid}.mfa'))

    # Load decoded states and calculate derivative
    posterior = []
    with open(f'out/{OGid}.tsv') as file:
        header = file.readline().split()
        for line in file:
            d = {key: float(value) for key, value in zip(header, line.split())}
            posterior.append(d['2'] + d['3'])
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
    ../../ortho_MSA/realign_hmmer/out/*
./decode.py
    ./out/*.tsv
"""