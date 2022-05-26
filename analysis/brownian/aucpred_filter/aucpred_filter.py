"""Remove sequences and regions from segments that do not pass quality filters."""

import re
import os

from src.utils import read_fasta


def load_posteriors(path):
    posteriors = {}
    with open(path) as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            ppid, p0, p1 = fields['ppid'], float(fields['0']), float(fields['1'])
            try:
                posteriors[ppid].append((p0, p1))
            except KeyError:
                posteriors[ppid] = [(p0, p1)]
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
spid_min = 20

# Load regions
OGid2regions = {}
with open('../aucpred_regions/out/regions.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGid, start, stop, disorder = fields['OGid'], int(fields['start']), int(fields['stop']), fields['disorder']
        try:
            OGid2regions[OGid].append((start, stop, disorder))
        except KeyError:
            OGid2regions[OGid] = [(start, stop, disorder)]

# Filter regions
record_sets = {i: [] for i in range(10, 35, 5)}
for OGid, regions in OGid2regions.items():
    msa = read_fasta(f'../../ortho_MSA/insertion_trim/out/{OGid}.afa')
    posteriors = load_posteriors(f'../../ortho_MSA/deletion_decode/out/{OGid}.tsv')

    for region in regions:
        # Get indices and length
        start, stop = region[0], region[1]
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
            if len(spids) >= spid_min and spid_filter(spids):
                record_sets[min_length].append((OGid, str(start), str(stop), disorder, ','.join(ppids)))

# Write records to file
if not os.path.exists('out/'):
    os.mkdir('out/')

for min_length, records in record_sets.items():
    with open(f'out/regions_{min_length}.tsv', 'w') as file:
        file.write('OGid\tstart\tstop\tdisorder\tppids\n')
        for record in records:
            file.write('\t'.join(record) + '\n')

"""
DEPENDENCIES
../../ortho_MSA/deletion_decode/decode.py
    ../../ortho_MSA/deletion_decode/out/*.tsv
../../ortho_MSA/insertion_trim/extract.py
    ../../ortho_MSA/insertion_trim/out/*.afa
../aucpred_regions/get_regions.py
    ../aucpred_regions/out/regions.tsv
"""