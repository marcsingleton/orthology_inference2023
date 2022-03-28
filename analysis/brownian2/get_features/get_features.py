"""Calculate features of segments in regions."""

import multiprocessing as mp
import os
import re
from collections import namedtuple

import features
from src.utils import read_fasta


Record = namedtuple('Record', ['OGid', 'start', 'stop', 'ppid', 'disorder', 'segment'])


def get_features(record):
    d = {'OGid': record.OGid, 'start': record.start, 'stop': record.stop, 'ppid': record.ppid}
    if not (len(record.segment) == 0 or 'X' in record.segment or 'U' in record.segment):
        d.update(features.get_features(record.segment))
    return d


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
ppid_regex = r'ppid=([A-Za-z0-9_]+)'

if __name__ == '__main__':
    # Load regions
    OGid2regions = {}
    with open('../aucpred_regions/out/regions.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            OGid, start, stop, disorder = line.split()
            try:
                OGid2regions[OGid].append((int(start), int(stop), disorder))
            except KeyError:
                OGid2regions[OGid] = [(int(start), int(stop), disorder)]

    # Extract segments
    args = []
    for OGid, regions in OGid2regions.items():
        msa = read_fasta(f'../insertion_trim/out/{OGid}.afa')
        msa = {re.search(ppid_regex, header).group(1): seq for header, seq in msa}

        for start, stop, disorder in regions:
            for ppid, seq in msa.items():
                segment = seq[start:stop].translate({ord('-'): None, ord('.'): None}).upper()
                args.append(Record(OGid, start, stop, ppid, disorder, segment))

    # Calculate features
    with mp.Pool(processes=num_processes) as pool:
        records = pool.map(get_features, args, chunksize=50)

    # Write features to file
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open('out/features.tsv', 'w') as file:
        if records:
            fields = list(records[0])
            file.write('\t'.join(fields) + '\n')
        for record in records:
            file.write('\t'.join(str(record.get(field, 'nan')) for field in fields) + '\n')

"""
DEPENDENCIES
../aucpred_regions/get_regions.py
    ../aucpred_regions/out/regions.tsv
../insertion_trim/extract.py
    ../insertion_trim/out/*.afa
"""