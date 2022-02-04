"""Create shuffled subsequences with identical amino acid composition and length distributions as the original classes."""

import os
import pandas as pd
import re
from itertools import accumulate
from random import shuffle

segment_dir = '../sample_segs/out/'  # Segment directory must end in /
type_name = 'ordered'  # Name of column denoting segment type

if not os.path.exists('out/'):
    os.mkdir('out/')

paths = filter(lambda x: re.match(r'segments_[0-9]+\.tsv', x), os.listdir(segment_dir))
for path in paths:
    # Load data and split into subsets
    segs = pd.read_csv(segment_dir + path, sep='\t', keep_default_na=False)
    T_seqs = segs.loc[segs[type_name], 'seq']
    F_seqs = segs.loc[~segs[type_name], 'seq']

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = path[j0+1:j1]

    # Create list of symbols and shuffle in place
    T_symlist = [sym for seq in T_seqs for sym in seq]
    F_symlist = [sym for seq in F_seqs for sym in seq]
    shuffle(T_symlist)
    shuffle(F_symlist)

    # Create splice indices
    T_lens = [len(seq) for seq in T_seqs]
    F_lens = [len(seq) for seq in F_seqs]
    T_accum = accumulate(T_lens)
    F_accum = accumulate(F_lens)

    # Slice shuffled lists
    shuffseq_T = [''.join(T_symlist[accum - size:accum]) for accum, size in zip(T_accum, T_lens)]
    shuffseq_F = [''.join(F_symlist[accum - size:accum]) for accum, size in zip(F_accum, F_lens)]

    # Convert to dataframe and save
    shuffseq = pd.DataFrame({type_name: segs[type_name], 'seq': shuffseq_T + shuffseq_F})
    shuffseq.to_csv(f'out/shuffseq_{i}.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../sample_segs/sample_seg.py
    ../sample_segs/out/segments_*.tsv
"""