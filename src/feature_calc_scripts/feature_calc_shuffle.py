"""Calculate features for different length cutoffs."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ, listdir
from sys import argv

dir = argv[1]  # Directory of segmented sequences must end in /
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    paths = filter(lambda x: x.endswith('.tsv'), listdir(dir))
    for path in paths:
        # Load data and subset
        df = pd.read_csv(dir + path, sep='\t', keep_default_na=False)
        df_con = df[df['conserved'] == True]
        df_div = df[df['conserved'] == False]

        # Get file index
        j0 = path.find('_')
        j1 = path.find('.tsv')
        i = path[j0 + 1:j1]

        # Compute features
        with mp.Pool(processes=num_processes) as pool:
            res_con = pd.DataFrame(pool.imap(seqfeat.feat_all, df_con['seq'], chunksize=50))
            res_div = pd.DataFrame(pool.imap(seqfeat.feat_all, df_div['seq'], chunksize=50))

        # Merge subsets and save
        res_all = pd.concat([res_con, res_div], keys=['con', 'div'])
        res_all.to_csv(f'features_{i}.tsv', sep='\t')
