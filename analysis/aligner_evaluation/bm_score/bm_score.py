"""Score alignments using SP and CS measures."""

import os
from itertools import combinations


def read_msf(path):
    with open(path) as file:
        header = False
        seqlines = {}
        for line in file:
            if line == '//\n':
                header = True
            elif header and line.rstrip('\n'):
                fields = line.split()
                seqid = fields[0]
                seqline = ''.join(fields[1:])
                try:
                    seqlines[seqid].append(seqline)
                except KeyError:
                    seqlines[seqid] = [seqline]
    seqs = {seqid: ''.join(seqlines) for seqid, seqlines in seqlines.items()}
    return seqs


def read_fasta(path):
    with open(path) as file:
        seqs = {}
        seqlines = []
        seqid = file.readline()[1:-1]
        line = file.readline()
        while line:
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seqs[seqid] = ''.join(seqlines)
            seqlines = []
            seqid = line[1:-1]
            line = file.readline()
    return seqs


def get_column_set(alignment, seqids=None, cutoff=1):
    if seqids is None:
        seqids = list(alignment)
    idx_seqs = [get_idx_seq(alignment[seqid]) for seqid in seqids]
    column_set = set()
    for col in set(zip(*idx_seqs)):
        if col.count(('-', None)) >= cutoff * len(col):
            continue
        column_set.add(col)
    return column_set


def get_pair_set(alignment, seqids=None, cutoff=1):
    if seqids is None:
        seqids = list(alignment)
    idx_seqs = [get_idx_seq(alignment[seqid]) for seqid in seqids]
    pair_set = set()
    for col in zip(*idx_seqs):
        if col.count(('-', None)) >= cutoff * len(col):
            continue
        for pair in combinations(filter(lambda x: x[1] != ('-', None), zip(seqids, col)), 2):
            pair_set.add(pair)
    return pair_set


def get_idx_seq(seq):
    i = 0
    idx_seq = []
    for sym in seq:
        if sym == '.' or sym == '-':
            idx_seq.append(('-', None))
        else:
            idx_seq.append((sym, i))
            i += 1
    return idx_seq


rows = []
aligners = [name for name in os.listdir('../bm_align/out/') if os.path.isdir('../bm_align/out/' + name)]
for aligner in aligners:
    for ref in ['11', '12', '20', '30', '40', '50']:
        file_ids = [name[:-4] for name in os.listdir(f'../bm_align/out/{aligner}/RV{ref}/') if name.endswith('.mfa')]
        for file_id in file_ids:
            ref_seq = read_msf(f'../../../data/bb3_release/RV{ref}/{file_id}.msf')
            test_seq = read_fasta(f'../bm_align/out/{aligner}/RV{ref}/{file_id}.mfa')
            seqids = list(ref_seq)

            ref_columns1 = get_column_set(ref_seq, seqids)
            ref_columns2 = get_column_set(ref_seq, seqids, 0.2)
            test_columns = get_column_set(test_seq, seqids)

            ref_pairs1 = get_pair_set(ref_seq, seqids)
            ref_pairs2 = get_pair_set(ref_seq, seqids, 0.2)
            test_pairs = get_pair_set(test_seq, seqids)

            rows.append({'aligner': aligner, 'ref': ref, 'file_id': file_id,
                         'ref_columns1': len(ref_columns1), 'test_columns1': len(test_columns & ref_columns1),
                         'ref_columns2': len(ref_columns2), 'test_columns2': len(test_columns & ref_columns2),
                         'ref_pairs1': len(ref_pairs1), 'test_pairs1': len(test_pairs & ref_pairs1),
                         'ref_pairs2': len(ref_pairs2), 'test_pairs2': len(test_pairs & ref_pairs2)})

if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/scores.tsv', 'w') as file:
    header = ['aligner', 'ref', 'file_id',
              'ref_columns1', 'test_columns1', 'ref_columns2', 'test_columns2',
              'ref_pairs1', 'test_pairs1', 'ref_pairs2', 'test_pairs2']
    file.write('\t'.join(header) + '\n')
    for row in rows:
        file.write('\t'.join([str(row[header]) for header in header]) + '\n')

"""
DEPENDENCIES
../../../data/bb3_release/
../bm_align/bm_align.py
    ../bm_align/bm_align/out/*
"""