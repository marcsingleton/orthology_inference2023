"""Calculate features for different length cutoffs."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ
from sys import argv

# Input variables
path = argv[1]  # Path to segmented sequences
key = argv[2]  # Key of column denoting subsequences class
T_idx = argv[3]  # Index of True class
F_idx = argv[4]  # Index of False class
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    df = pd.read_csv(path, sep='\t', keep_default_na=False)

    for i in [2 ** x for x in range(6)]:
        # Load data and subset
        df_cutoff = df[df['seq'].map(lambda x: len(x) >= i)]  # Select entries where sequence is above length threshold
        df_T = df_cutoff[df_cutoff[key] == True].sample(8500)
        df_F = df_cutoff[df_cutoff[key] == False].sample(8500)
        df_seq = pd.concat([df_T, df_F])
        df_seq.to_csv(f'sequences_{i}.tsv', sep='\t')

        # Compute features
        with mp.Pool(processes=num_processes) as pool:
            res_T = pd.DataFrame(pool.imap(seqfeat.feat_all, df_T['seq'], chunksize=50))
            res_F = pd.DataFrame(pool.imap(seqfeat.feat_all, df_F['seq'], chunksize=50))

        # Merge subsets and save
        res_all = pd.concat([res_T, res_F], keys=[T_idx, F_idx])
        res_all.to_csv(f'features_{i}.tsv', sep='\t')
