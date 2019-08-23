"""Create shuffled conserved and diverged subsequence with identical amino acid composition and length distirbutions."""

import os
import pandas as pd
import re
from itertools import accumulate
from random import shuffle

sequence_dir = '../feature_calc/'
paths = filter(lambda x: re.match('sequences_[0-9]+\.tsv', x), os.listdir(sequence_dir))
for path in paths:
    # Load data and split into subsets
    df = pd.read_csv(sequence_dir + path, sep='\t', keep_default_na=False, index_col=0)
    seq_con = df[df['conserved'] == True]['seq']
    seq_div = df[df['conserved'] == False]['seq']

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = path[j0+1:j1]

    # Create list of symbols and shuffle in place
    shufflist_con = [sym for seq in seq_con for sym in seq]
    shufflist_div = [sym for seq in seq_div for sym in seq]
    shuffle(shufflist_con)
    shuffle(shufflist_div)

    # Create splice indices
    lens_con = [len(seq) for seq in seq_con]
    lens_div = [len(seq) for seq in seq_div]
    accum_con = accumulate(lens_con)
    accum_div = accumulate(lens_div)

    # Slice shuffled lists
    shuffseq_con = [''.join(shufflist_con[accum-size:accum]) for accum, size in zip(accum_con, lens_con)]
    shuffseq_div = [''.join(shufflist_div[accum-size:accum]) for accum, size in zip(accum_div, lens_div)]

    # Convert to dataframe and save
    df_shuff = pd.DataFrame({'conserved': [True] * len(shuffseq_con) + [False] * len(shuffseq_div),
                             'seq': shuffseq_con + shuffseq_div})
    df_shuff.to_csv(f'shuffseq_{i}.tsv', sep='\t')

"""
DEPENDENCIES
../feature_calc/feature_calc.py
    ../feature_calc/sequences_?.tsv
"""