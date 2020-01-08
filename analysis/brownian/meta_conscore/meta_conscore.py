"""Calculate conservation scores for IDR segments."""

import pandas as pd
from Bio.SubsMat import MatrixInfo
from itertools import combinations

# Input variables
path = '../segment_avg/segment_avg.tsv'
matrix = MatrixInfo.blosum50

# Initialize
scores_blosum = {}
scores_match = {}
segs = pd.read_csv(path, sep='\t', keep_default_na=False)

for key, block in segs.groupby('block_id'):
    MSA = block['seq'].array
    rows, cols = len(MSA), len(MSA[0])
    score_blosum = 0
    score_match = 0
    for j in range(cols):
        col = [MSA[i][j] for i in range(rows)]  # Create column list of symbols
        for pair in combinations(col, 2):
            score_blosum += matrix.get(pair, -2)
            score_match += 1 if pair[0] != '-' and pair[0] == pair[1] else 0  # Do not count gap-gap matches

    scores_blosum[key] = score_blosum / (rows * (rows - 1) * cols / 2)
    scores_match[key] = score_match / (rows * (rows - 1) * cols / 2)

# Save block_id: score dictionaries as tsv
with open('scores_blosum.tsv', 'w') as file:
    file.write('block_id\tblosum_score\n')
    for block_id, blosum_score in scores_blosum.items():
        file.write(f'{block_id}\t{blosum_score}\n')
with open('scores_match.tsv', 'w') as file:
    file.write('block_id\tmatch_score\n')
    for block_id, match_score in scores_match.items():
        file.write(f'{block_id}\t{match_score}\n')

"""
DEPENDENCIES
../segment_avg/segment_avg.py
    ../segment_avg/segment_avg.tsv
"""