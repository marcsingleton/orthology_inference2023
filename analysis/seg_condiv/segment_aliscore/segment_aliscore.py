"""Segment alignments by alignment score and write subsequences to dataframe."""

import gzip
import pandas as pd
from Bio import AlignIO
from Bio.SubsMat import MatrixInfo
from itertools import combinations
from scipy import ndimage

path = '../../EggNOGv5_validation/filter_count/out/7214_noX_members/10_10_members.tsv'
dir = '../../EggNOGv5_validation/filter_unknown_realign/out/align/'
matrix = MatrixInfo.blosum50
thresh = 0

segs = []  # Segment data with raw sequences
seg_num = 0  # Counter for numbering rows

with open(path) as file:
    for line in file:
        fields = line.split('\t')
        ali_id = fields[1]
        with gzip.open(dir + ali_id + '.raw_alg.faa.gz', 'rt') as file_MSA:  # read, text mode
            MSA = AlignIO.read(file_MSA, 'fasta')
            row_len = len(MSA[0, :])
            column_len = len(MSA)

            # Calculate score of each column
            column_scores = []
            for i in range(row_len):
                score = 0
                for pair in combinations(MSA[:, i], 2):
                    score += matrix.get(pair, -2)
                column_scores.append(score)
            column_scores = ndimage.gaussian_filter1d(column_scores, 2)

            # Find bounds
            bounds = []
            conserved_prev = column_scores[0] > thresh
            bound_start = 0
            for i, score in enumerate(column_scores):
                conserved_curr = score > thresh
                if conserved_curr is not conserved_prev:
                    bounds.append(((bound_start, i), conserved_prev))
                    bound_start = i
                conserved_prev = conserved_curr
            bounds.append(((bound_start, i + 1), conserved_curr))  # Bound for final segment

            # Create dataframe rows
            for bound, conserved in bounds:
                for record in MSA:
                    seq_raw = str(record.seq[slice(*bound)])
                    seq_ungap = seq_raw.translate({ord('-'): None})
                    segs.append({'ali_id': ali_id, 'seq_id': record.id, 'seg_id': hex(seg_num)[2:].zfill(8),
                                 'bound': bound, 'conserved': conserved, 'seq': seq_raw})
                    seg_num += 1

df = pd.DataFrame(segs)
df.to_csv('segment_aliscore.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../../EggNOGv5_validation/filter_count/filter_count.py
    ../../EggNOGv5_validation/filter_count/out/7214_noX_members/10_10_members.tsv
../../EggNOGv5_validation/filter_unknown_realign/filter_unknown_realign.py
    ../../EggNOGv5_validation/filter_unknown_realign/out/align/
"""