"""Calculate features for different length cutoffs."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ
from sys import argv

# Input variables
path = argv[1]  # Path to segmented sequences
key = argv[2]  # Key of column denoting subsequences class
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    segs = pd.read_csv(path, sep='\t', keep_default_na=False)

    for i in [2 ** x for x in range(6)]:
        # Load data and subset
        segs_lt = segs[segs['seq'].map(lambda x: len(x) >= i)]  # Select entries where sequence is above length threshold
        sample_T = segs_lt[segs_lt[key] == True].sample(8500)
        sample_F = segs_lt[segs_lt[key] == False].sample(8500)
        sample = pd.concat([sample_T, sample_F])
        sample.to_csv(f'sequences_{i}.tsv', sep='\t')

        # Compute features
        with mp.Pool(processes=num_processes) as pool:
            features_T = pd.DataFrame(pool.imap(seqfeat.feat_all, sample_T['seq'], chunksize=50))
            features_F = pd.DataFrame(pool.imap(seqfeat.feat_all, sample_F['seq'], chunksize=50))

        # Merge subsets and save
        features = pd.concat([features_T, features_F])
        features.set_index(sample[key], inplace=True)
        features.to_csv(f'features_{i}.tsv', sep='\t')
