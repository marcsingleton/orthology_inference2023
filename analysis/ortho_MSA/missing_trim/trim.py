"""Segment state 2 posterior decodings into contiguous intervals."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
from src.utils import read_fasta

ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
posterior_high = 0.95
posterior_low = 0.5

# Load OGids
OGids = []
with open('../realign_hmmer/out/errors.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGid, error_flag = fields['OGid'], fields['error_flag']
        if error_flag == 'False':
            OGids.append(OGid)

if not os.path.exists('out/trims/'):
    os.mkdir('out/trims/')

rows = []
for OGid in OGids:
    msa = {}
    for header, seq in read_fasta(f'../insertion_trim/out/trims/{OGid}.afa'):
        ppid = re.search(ppid_regex, header).group(1)
        msa[ppid] = seq

    posteriors = []
    with open(f'out/posteriors/{OGid}.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        fields = {key: value for key, value in zip(field_names, file.readline().rstrip('\n').split('\t'))}
        ppid0, p2 = fields['ppid'], float(fields['2'])
        posterior = [p2]
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            ppid, p2 = fields['ppid'], float(fields['2'])
            if ppid == ppid0:
                posterior.append(p2)
            else:
                posteriors.append((ppid0, np.array(posterior)))
                ppid0, posterior = ppid, [p2]
        posteriors.append((ppid0, np.array(posterior)))  # Append last posterior in file

    trims = []
    for ppid, posterior in posteriors:
        slices = []
        for s, in ndimage.find_objects(ndimage.label(posterior >= posterior_high)[0]):
            seq = msa[ppid]

            start = s.start
            while (start-1 >= 0) and (posterior[start-1] >= posterior_low) and (seq[start-1] not in ['-', '.']):
                start -= 1
            stop = s.stop
            while (stop+1 < len(seq)) and (posterior[stop+1] >= posterior_low) and (seq[stop+1] not in ['-', '.']):
                stop += 1

            slices.append((start, stop+1))
        trims.append((ppid, slices))

    with open(f'out/trims/{OGid}.tsv', 'w') as file:
        file.write('ppid\tslices\n')
        for ppid, slices in trims:
            slicestring = ','.join([f'{start}-{stop}' for start, stop in slices])
            file.write(f'{ppid}\t{slicestring}\n')

    for ppid, slices in trims:
        for start, stop in slices:
            rows.append({'OGid': OGid, 'colnum': len(posteriors[0][1]), 'ppid': ppid, 'start': start, 'stop': stop})

# Plot stats
df = pd.DataFrame(rows)
df.to_csv('out/trim_stats.tsv', sep='\t', index=False)

df['length'] = df['stop'] - df['start']
df['length_ratio'] = df['length'] / df['colnum']
ppid_groups = df.groupby(['OGid', 'ppid'])  # Use combined key in case ppid in multiple OGids
OGid_groups = df.groupby(['OGid'])

# Number of OGs with trims
values = [len(set(OGids)) - df['OGid'].nunique(), df['OGid'].nunique()]
labels = [f'{label}\n{value:,}' for label, value in zip(['w/o trims', 'w/ trims'], values)]
fig, ax = plt.subplots()
ax.pie(values, labels=labels, labeldistance=1.15)
ax.set_title('OGs w/ and w/o trims')
fig.savefig('out/pie_trims.png')
plt.close()

# Distribution of number sequences with trims in OGs
counts = OGid_groups['ppid'].nunique().value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of unique sequences with trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-seqnum.png')
plt.close()

# Distribution of number of trims in sequences
counts = ppid_groups.size().value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of trims in sequence')
ax.set_ylabel('Number of sequences')
fig.savefig('out/hist_OGnum-seqnum.png')
plt.close()

# Distribution of length of trims
fig, ax = plt.subplots()
ax.hist(df['length'], bins=100)
ax.set_xlabel('Length of trim')
ax.set_ylabel('Number of trims')
fig.savefig('out/hist_trimnum-length.png')
ax.set_yscale('log')
fig.savefig('out/hist_trimnum-length_log.png')
plt.close()

# Distribution of length ratio of trims
fig, ax = plt.subplots()
ax.hist(df['length_ratio'], bins=50)
ax.set_xlabel('Length ratio of trim')
ax.set_ylabel('Number of trims')
fig.savefig('out/hist_trimnum-ratio.png')
plt.close()

"""
NOTES
While unlikely, it is possible for a protein to be in multiple OGids at this point since during gene clustering, some
proteins were in multiple gene groups.

DEPENDENCIES
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/errors.tsv
./decode.py/
    ./out/posteriors/*.tsv
"""