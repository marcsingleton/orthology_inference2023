"""Segment trimmed alignments into regions by averaging raw AUCpreD scores."""

import os
import re

import numpy as np
from src.utils import read_fasta


def load_scores(path):
    with open(path) as file:
        scores = []
        for line in file:
            if not line.startswith('#'):
                score = line.split()[3]
                scores.append(score)
    return scores


def gaussian_filter(array, sigma):
    # Make stack of Gaussian kernels
    radius = int(4 * sigma + 0.5)  # Truncate filter at 4 standard deviations rounded to nearest integer
    x = np.stack([np.arange(-radius, radius+1) for _ in range(array.shape[0])])
    kernel = np.exp(-x**2 / (2 * sigma**2))

    # Apply filter ignoring masked values
    padded = np.ma.masked_invalid(np.pad(array, [(0, 0), (radius, radius)], mode='edge'))
    mean_array = np.zeros(array.shape[1])
    var_array = np.zeros(array.shape[1])
    for j in range(array.shape[1]):  # Output has as many columns as input; slicing always grabs correct window even though actual centers are offset
        window = padded[:, j:j+2*radius+1]
        weight = (~window.mask * kernel).sum()
        mean = (window * kernel).sum() / weight
        var = ((window - mean) ** 2 * kernel).sum() / weight

        mean_array[j] = mean
        var_array[j] = var

    return mean_array, var_array


threshold = 0.5
ppid_regex = r'ppid=([A-Za-z0-9_]+)'

records = []
for OGid in os.listdir('out/raw/'):
    # Load MSA
    msa = read_fasta(f'../insertion_trim/out/{OGid}.mfa')
    msa = {re.search(ppid_regex, header).group(1): seq for header, seq in msa}

    # Map outputs to MSA columns
    ppids = {path.split('.')[0] for path in os.listdir(f'out/raw/{OGid}/')}
    if set(msa) != ppids:
        print(f'{OGid} has fewer predictions than sequences. Skipping segmentation.')
        continue

    mapped = np.full((len(ppids), max([len(seq) for seq in msa.values()])), np.nan)
    for i, ppid in enumerate(ppids):
        scores = load_scores(f'out/raw/{OGid}/{ppid}.diso_noprof')
        idx = 0
        for j, sym in enumerate(msa[ppid]):
            if sym not in ['-', '.']:
                mapped[i, j] = scores[idx]
                idx += 1
    mapped = np.ma.masked_invalid(mapped)

    # Extract regions
    mean, var = gaussian_filter(mapped, 2)
    binary = mean >= threshold
    regions, value0, idx0 = [], binary[0], 0
    for idx, value in enumerate(binary):
        if value != value0:
            regions.append((idx0, idx, value0))
            value0, idx0 = value, idx
    regions.append((idx0, idx+1, value0))

    # Write regions as records
    for start, stop, disorder in regions:
        records.append((OGid, str(start), str(stop), str(disorder)))

# Write segments to file
with open('out/regions.tsv', 'w') as file:
    fields = ['OGid', 'start', 'stop', 'disorder']
    file.write('\t'.join(fields) + '\n')
    for record in records:
        file.write('\t'.join(record) + '\n')

"""
OUTPUT
3423 has fewer predictions than sequences. Skipping segmentation.

DEPENDENCIES
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
./run_aucpred.py
    ./out/raw/*/*.diso_noprof
"""