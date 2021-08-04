"""Segment trimmed alignments by averaging raw AUCpreD scores."""

import os
import re

import numpy as np


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


def load_scores(path):
    with open(path) as file:
        scores = []
        for line in file:
            if not line.startswith('#'):
                score = line.split()[3]
                scores.append(score)
    return scores


THRESHOLD = 0.5
PPID_REGEX = r'ppid=([A-Za-z0-9_]+)'

records = []
for OGid in os.listdir('out/raw/'):
    # Load MSA
    msa = load_msa(f'../trim_extract/out/{OGid}.mfa')
    msa = {re.search(PPID_REGEX, header).group(1): seq for header, seq in msa}

    # Map outputs to MSA columns
    ppids = set([path.split('.')[0] for path in os.listdir(f'out/raw/{OGid}/')])
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

    # Extract regions
    binary = np.nanmean(mapped, axis=0) >= THRESHOLD
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
with open('out/segments.tsv', 'w') as file:
    fields = ['OGid', 'start', 'stop', 'disorder']
    file.write('\t'.join(fields) + '\n')
    for record in records:
        file.write('\t'.join(record) + '\n')

"""
DEPENDENCIES
../trim_extract/trim_extract.py
    ../trim_extract/out/*.mfa
./run.py
    ./out/raw/*/*.diso_noprof
"""