"""Trim alignments, removing spurious sequences and regions."""

import os

import numpy as np
import skbio
import tensorflow as tf
from src.utils import read_fasta
from src.brownian2.trim import get_slices

alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', '-']
sym2idx = {sym: i for i, sym in enumerate(alphabet)}
ppid_regex = r'ppid=([A-Za-z0-9_.]+)'

posterior_high = 0.9
posterior_low = 0.005
gradient_high = np.inf
gradient_low = np.inf

OGids = []
with open('../OG_filter/out/OG_filter.tsv') as file:
    header = file.readline().rstrip('\n').split('\t')
    field2idx = {field: i for i, field in enumerate(header)}
    for line in file:
        fields = line.rstrip('\n').split('\t')
        OGid = fields[field2idx['OGid']]
        OGids.append(OGid)

model = tf.keras.models.load_model('../cnn_fit/out/model.h5')

if not os.path.exists('out/'):
    os.mkdir('out')

for OGid in OGids:
    try:
        msa1 = read_fasta(f'../align_fastas1/out/{OGid}.afa')
    except FileNotFoundError:
        msa1 = read_fasta(f'../align_fastas2/out/{OGid}.afa')

    # Convert alignment to indices and make profile
    msa2 = []
    for i, (header, seq1) in enumerate(msa1):
        seq2 = [sym2idx.get(sym, sym2idx['X']) for sym in seq1]  # All non-standard symbols mapped to X
        msa2.append(seq2)
    msa2 = tf.keras.utils.to_categorical(msa2, len(alphabet))
    profile = msa2.sum(axis=0) / len(msa2)

    # Identify trims
    trims_array = np.full((len(msa1), len(msa1[0][1])), False)
    for i, seq in enumerate(msa2):
        output = tf.squeeze(model([np.expand_dims(profile, 0), np.expand_dims(seq, 0)]))  # Expand and contract dims
        gradient = np.gradient(output)
        slices = get_slices(msa1, output, gradient, posterior_high, posterior_low, gradient_high, gradient_low)
        for s in slices:
            trims_array[i, s] = True

    # Identify gaps
    gaps_array = np.full((len(msa1), len(msa1[0][1])), False)
    for i, (_, seq) in enumerate(msa1):
        for j, sym in enumerate(seq):
            if sym == '-':
                gaps_array[i, j] = True

    # Identify consensus columns
    scores = np.logical_or(trims_array, gaps_array).sum(axis=0)
    rf = ['.' if score == len(msa1) else 'x' for score in scores]  # Metadata for marking consensus columns in profile HMM

    # Write to file
    msa3 = skbio.TabularMSA([skbio.Protein(seq, metadata={'description': header}) for header, seq in msa1],
                            positional_metadata={'RF': rf})
    msa3.write(f'out/{OGid}.sto', 'stockholm')

"""
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.afa
../align_fastas2/align_fastas2.py
    ../align_fastas2/out/*.afa
../cnn_fit/cnn_fit.py
    ../cnn_fit/out/model.h5
../OG_filter/OG_filter.py
    ../OG_filter/out/OG_filter.tsv
"""