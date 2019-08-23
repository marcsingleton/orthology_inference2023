"""Create shuffled subsequences with identical amino acid composition and length distirbutions as the original classes."""

import os
import pandas as pd
import re
from itertools import accumulate
from random import shuffle
from sys import argv

sequence_dir = '../feature_calc/'
key = argv[1]  # Key of column denoting subsequences class
paths = filter(lambda x: re.match('sequences_[0-9]+\.tsv', x), os.listdir(sequence_dir))
for path in paths:
    # Load data and split into subsets
    df = pd.read_csv(sequence_dir + path, sep='\t', keep_default_na=False, index_col=0)
    seq_T = df[df[key] == True]['seq']
    seq_F = df[df[key] == False]['seq']

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = path[j0+1:j1]

    # Create list of symbols and shuffle in place
    shufflist_T = [sym for seq in seq_T for sym in seq]
    shufflist_F = [sym for seq in seq_F for sym in seq]
    shuffle(shufflist_T)
    shuffle(shufflist_F)

    # Create splice indices
    lens_T = [len(seq) for seq in seq_T]
    lens_F = [len(seq) for seq in seq_F]
    accum_T = accumulate(lens_T)
    accum_F = accumulate(lens_F)

    # Slice shuffled lists
    shuffseq_T = [''.join(shufflist_T[accum - size:accum]) for accum, size in zip(accum_T, lens_T)]
    shuffseq_F = [''.join(shufflist_F[accum - size:accum]) for accum, size in zip(accum_F, lens_F)]

    # Convert to dataframe and save
    df_shuff = pd.DataFrame({key: [True] * len(shuffseq_T) + [False] * len(shuffseq_F),
                             'seq': shuffseq_T + shuffseq_F})
    df_shuff.to_csv(f'shuffseq_{i}.tsv', sep='\t')
