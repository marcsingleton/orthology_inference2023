"""Generate samples from ASR distributions."""

import os
import random

import numpy as np

random.seed(930715)  # Set seed to make results consistent
syms = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
idx2sym = {i: sym for i, sym in enumerate(syms)}

if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = [path[:4] for path in os.listdir('../asr_root/out/') if path.endswith('_aa.npy')]
for OGid in sorted(OGids):  # Ensure consistent order and thus consistent samples
    # Generate amino acid sequence
    aa_dist = np.load(f'../asr_root/out/{OGid}_aa.npy')
    aa_marginal = aa_dist.sum(axis=0)

    seqs = []
    for _ in range(1000):
        seq = []
        for i in range(aa_marginal.shape[1]):
            idx, = random.choices(list(range(20)), aa_marginal[:, i])
            sym = idx2sym[idx]
            seq.append(sym)
        seqs.append(seq)

    # Generate indels (if applicable)
    if os.path.exists(f'../asr_root/out/{OGid}_indel.npy'):
        indel_dist = np.load(f'../asr_root/out/{OGid}_indel.npy')
        indel_marginal = indel_dist.sum(axis=0)

        idx2indel = {}
        with open(f'../asr_indel/out/{OGid}.tsv') as file:
            file.readline()  # Skip header
            for line in file:
                idx, start, stop = [int(field) for field in line.split()]
                idx2indel[idx] = (start, stop)

        for seq in seqs:
            for i in range(indel_marginal.shape[1]):
                p = random.random()
                if p > indel_marginal[0, i]:
                    start, stop = idx2indel[i]
                    seq[start:stop] = (stop-start)*['-']

    with open(f'out/{OGid}_sample.mfa', 'w') as file:
        for idx, seq in enumerate(seqs):
            seq = ''.join(seq)
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]) + '\n'
            file.write(f'>seq{idx}\n' + seqstring)

"""
NOTES
The output of the ancestral reconstructions gives a probability for each rate and symbol combination. One method of
sampling sequences and rates is directly from this posterior distribution. Thus, every draw for a position yields both a
symbol and rate, where the rates are the categories fit by the model. While this method would capture the uncertainty in
the rate inference, I feel it's both overly complex and a poor interpretation of the model in this context. The gamma
distribution actually describes continuous variation in rates even though it's discretized for computational
convenience here. It's also more realistic to expect the biological variation in rates to be continuous. Thus, rather
than sampling the rate of site as "invariant 10% of the time, slow 25% of the time...", it's more appropriate to take
the rate as a posterior average over rate categories. I feel this will yield better simulated alignments since the
rate variation will manifest as a continuum rather than a single site being fast in some simulations and slow in others.

DEPENDENCIES
../asr_root/aa.py
    ../asr_root/out/*_aa.npy
../asr_root/indel.py
    ../asr_root/out/*_indel.npy
"""