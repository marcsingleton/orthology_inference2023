"""Calculate features of set of shuffled sequences."""

import multiprocessing as mp
import os
import pandas as pd
import src.seg_scripts.seqfeat as seqfeat

# Input variables
segment_dir = '../segment_shuffle/out/'  # Directory of segmented sequences must end in /
type_name = 'ordered'  # Name of column denoting segment type
num_processes = int(os.environ['SLURM_NTASKS'])

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    paths = filter(lambda x: x.endswith('.tsv'), os.listdir(segment_dir))
    for path in paths:
        # Load data and subset
        segs = pd.read_csv(segment_dir + path, sep='\t', keep_default_na=False)
        T_segs = segs[segs[type_name]]
        F_segs = segs[~segs[type_name]]  # ~ is bitwise NOT operator; it interacts properly with numpy objects but not Python booleans

        # Get file index
        j0 = path.find('_')
        j1 = path.find('.tsv')
        i = path[j0 + 1:j1]

        # Compute features
        with mp.Pool(processes=num_processes) as pool:
            T_features = pd.DataFrame(pool.imap(seqfeat.feat_all, T_segs['seq'], chunksize=50))
            F_features = pd.DataFrame(pool.imap(seqfeat.feat_all, F_segs['seq'], chunksize=50))

        # Merge subsets and save
        features = pd.concat([T_features, F_features])
        features.set_index(segs[type_name], inplace=True)
        features.to_csv(f'out/features_{i}.tsv', sep='\t')

"""
DEPENDENCIES
../../../src/seg_scripts/seqfeat.py
../segment_shuffle/segment_shuffle.py
    ../segment_shuffle/out/shuffseq_*.tsv
"""