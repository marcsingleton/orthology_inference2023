"""Calculate features for different length cutoffs."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ, listdir
from sys import argv

# Input variables
dir = argv[1]  # Directory of segmented sequences must end in /
key = argv[2]  # Key of column denoting subsequences class
T_idx = argv[3]  # Index of True class in sentence case
F_idx = argv[4]  # Index of False class in sentence case
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    paths = filter(lambda x: x.endswith('.tsv'), listdir(dir))
    for path in paths:
        # Load data and subset
        df = pd.read_csv(dir + path, sep='\t', keep_default_na=False)
        df_T = df[df[key] == True]
        df_F = df[df[key] == False]

        # Get file index
        j0 = path.find('_')
        j1 = path.find('.tsv')
        i = path[j0 + 1:j1]

        # Compute features
        with mp.Pool(processes=num_processes) as pool:
            res_T = pd.DataFrame(pool.imap(seqfeat.feat_all, df_T['seq'], chunksize=50))
            res_F = pd.DataFrame(pool.imap(seqfeat.feat_all, df_F['seq'], chunksize=50))

        # Merge subsets and save
        res_all = pd.concat([res_T, res_F], keys=[T_idx, F_idx])
        res_all.to_csv(f'features_{i}.tsv', sep='\t')
