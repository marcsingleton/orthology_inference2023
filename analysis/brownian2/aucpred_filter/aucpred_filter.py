"""Remove sequences and regions from segments that do not pass quality filters."""

import re
import os

from src.utils import read_fasta


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


threshold = 0.99
ppid_regex = r'ppid=([A-Za-z0-9_]+)'
spid_regex = r'spid=([a-z]+)'
alphabet = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}

# Load regions
OGid2regions = {}
with open('../aucpred_regions/out/regions.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder = line.split()
        try:
            OGid2regions[OGid].append((start, stop, disorder))
        except KeyError:
            OGid2regions[OGid] = [(start, stop, disorder)]

# Filter regions
record_sets = {i: [] for i in range(10, 35, 5)}
for OGid, regions in OGid2regions.items():
    msa = read_fasta(f'../insertion_trim/out/{OGid}.mfa')
    posteriors = load_posteriors(f'../deletion_decode/out/{OGid}.tsv')

    for region in regions:
        # Get indices and length
        start, stop = int(region[0]), int(region[1])
        length = stop - start
        disorder = region[2]

        # Extract and filter segments
        segment_sets = {i: [] for i in record_sets}
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
                if sym not in alphabet:
                    is_standard = False
                    break
            posterior = posteriors[ppid][start:stop]
            for min_length, segments in segment_sets.items():
                if length >= min_length and is_standard and all([p1 < threshold for _, p1 in posterior]):
                    segments.append((ppid, spid))

        # Filter by phylogenetic diversity
        for min_length, segments in segment_sets.items():
            ppids = [ppid for ppid, _ in segments]
            spids = {spid for _, spid in segments}
            if len(spids) >= 20 and spid_filter(spids):
                record_sets[min_length].append((OGid, str(start), str(stop), str(disorder), ','.join(ppids)))

# Write records to file
if not os.path.exists('out/'):
    os.mkdir('out/')

for min_length, records in record_sets.items():
    with open(f'out/regions_{min_length}.tsv', 'w') as file:
        fields = ['OGid', 'start', 'stop', 'disorder', 'ppids']
        file.write('\t'.join(fields) + '\n')
        for record in records:
            file.write('\t'.join(record) + '\n')

"""
DEPENDENCIES
../aucpred_regions/get_regions.py
    ../aucpred_regions/out/regions.tsv
../deletion_decode/decode.py
    ../deletion_decode/out/*.tsv
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
"""