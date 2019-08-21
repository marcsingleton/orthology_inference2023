"""Calculate features for different length cutoffs."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ
from sys import argv

path = argv[1]  # Path to segmented sequences
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    df = pd.read_csv(path, sep='\t', keep_default_na=False)

    for i in [2 ** x for x in range(6)]:
        # Load data and subset
        df_cutoff = df[df['seq'].map(lambda x: len(x) >= i)]  # Select entries where sequence is above length threshold
        df_con = df_cutoff[df_cutoff['conserved'] == True].sample(8500)
        df_div = df_cutoff[df_cutoff['conserved'] == False].sample(8500)
        df_seq = pd.concat([df_con, df_div])
        df_seq.to_csv(f'sequences_{i}.tsv', sep='\t')

        # Compute features
        with mp.Pool(processes=num_processes) as pool:
            res_con = pd.DataFrame(pool.imap(seqfeat.feat_all, df_con['seq'], chunksize=50))
            res_div = pd.DataFrame(pool.imap(seqfeat.feat_all, df_div['seq'], chunksize=50))

        # Merge subsets and save
        res_all = pd.concat([res_con, res_div], keys=['con', 'div'])
        res_all.to_csv(f'features_{i}.tsv', sep='\t')
