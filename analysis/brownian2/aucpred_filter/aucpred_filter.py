"""Remove sequences and regions from segments that do not pass quality filters."""

import re
import os

import scipy.stats as stats


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


def load_posteriors(path):
    posteriors = {}
    with open(path) as file:
        file.readline()  # Skip header
        for line in file:
            ppid, p0, p1 = line.split()
            try:
                posteriors[ppid].append((float(p0), float(p1)))
            except KeyError:
                posteriors[ppid] = [(float(p0), float(p1))]
    return posteriors


def spid_filter(spids):
    conditions = [({'dnov', 'dvir'}, 1),
                  ({'dmoj', 'dnav'}, 1),
                  ({'dinn', 'dgri', 'dhyd'}, 2),
                  ({'dgua', 'dsob'}, 1),
                  ({'dbip', 'dana'}, 1),
                  ({'dser', 'dkik'}, 1),
                  ({'dele', 'dfik'}, 1),
                  ({'dtak', 'dbia'}, 1),
                  ({'dsuz', 'dspu'}, 1),
                  ({'dsan', 'dyak'}, 1),
                  ({'dmel'}, 1),
                  ({'dmau', 'dsim', 'dsec'}, 1)]
    return all([len(spids & group) >= num for group, num in conditions])


min_length = 30
threshold = 0.99
ppid_regex = r'ppid=([A-Za-z0-9_]+)'
spid_regex = r'spid=([a-z]+)'

# Load regions
OGid2regions = {}
with open('../aucpred_segment/out/segments.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder = line.split()
        try:
            OGid2regions[OGid].append((start, stop, disorder))
        except KeyError:
            OGid2regions[OGid] = [(start, stop, disorder)]

# Filter regions
records = []
for OGid, regions in OGid2regions.items():
    msa = load_msa(f'../insertion_trim/out/{OGid}.mfa')
    posteriors = load_posteriors(f'../deletion_decode/out/{OGid}.tsv')

    for region in regions:
        # Get indices and length
        start, stop = int(region[0]), int(region[1])
        length = stop - start
        disorder = region[2]

        # Extract and filter segments
        segments = []
        for header, seq in msa:
            # Extract
            ppid = re.search(ppid_regex, header).group(1)
            spid = re.search(spid_regex, header).group(1)
            segment = seq[start:stop].upper()

            # Filter
            length, is_standard = 0, True
            for sym in segment:
                if sym not in ['-', '.']:
                    length += 1
                if sym in ['X', 'U']:
                    is_standard = False
                    break
            posterior = posteriors[ppid][start:stop]
            if length >= min_length and is_standard and all([p1 < threshold for _, p1 in posterior]):
                segments.append((ppid, spid))

        # Filter by phylogenetic diversity
        ppids = [ppid for ppid, _ in segments]
        spids = set([spid for _, spid in segments])
        if len(spids) >= 20 and spid_filter(spids):
            records.append((OGid, str(start), str(stop), str(disorder), ','.join(ppids)))

# Write records to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/segments.tsv', 'w') as file:
    fields = ['OGid', 'start', 'stop', 'disorder', 'ppids']
    file.write('\t'.join(fields) + '\n')
    for record in records:
        file.write('\t'.join(record) + '\n')

"""
DEPENDENCIES
../aucpred_segment/segment.py
    ../aucpred_segment/out/segments.tsv
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
"""