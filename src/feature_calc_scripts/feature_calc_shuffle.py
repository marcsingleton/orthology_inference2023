"""Calculate features for different length cutoffs."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ, listdir
from sys import argv

# Input variables
dir = argv[1]  # Directory of segmented sequences must end in /
key = argv[2]  # Key of column denoting subsequences class
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    paths = filter(lambda x: x.endswith('.tsv'), listdir(dir))
    for path in paths:
        # Load data and subset
        segs = pd.read_csv(dir + path, sep='\t', keep_default_na=False)
        segs_T = segs[segs[key] == True]
        segs_F = segs[segs[key] == False]

        # Get file index
        j0 = path.find('_')
        j1 = path.find('.tsv')
        i = path[j0 + 1:j1]

        # Compute features
        with mp.Pool(processes=num_processes) as pool:
            features_T = pd.DataFrame(pool.imap(seqfeat.feat_all, segs_T['seq'], chunksize=50))
            features_F = pd.DataFrame(pool.imap(seqfeat.feat_all, segs_F['seq'], chunksize=50))

        # Merge subsets and save
        features = pd.concat([features_T, features_F])
        features.set_index(segs[key], inplace=True)
        features.to_csv(f'features_{i}.tsv', sep='\t')
