"""Calculate features of set of shuffled sequences."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ, listdir
from sys import argv

# Input variables
segment_dir = argv[1]  # Directory of segmented sequences must end in /
type_name = argv[2]  # Name of column denoting segment type
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    paths = filter(lambda x: x.endswith('.tsv'), listdir(segment_dir))
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
        features.to_csv(f'features_{i}.tsv', sep='\t')
