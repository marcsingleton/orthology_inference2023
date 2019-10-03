"""Calculate features for different length cutoffs while maintaining block identity."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ

# Input variables
path = '../segment_avg/segment_avg.tsv'  # Path to segmented sequences
key = 'ordered'  # Key of column denoting block class
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    df = pd.read_csv(path, sep='\t', keep_default_na=False)

    for i in [2 ** x for x in range(6)]:
        # Load data and subset
        blocks_seq = []
        blocks_id = []
        blocks_class = []
        blocks_species = []
        for _, block in df.groupby('block_id'):
            seqs = block['seq'].map(lambda x: x.translate({ord('-'): None}))
            if seqs.map(lambda x: len(x) >= i).all():  # All subsequences in block must meet length cutoff
                blocks_seq.extend(seqs)
                blocks_id.extend(block['block_id'])
                blocks_class.extend(block[key])
                blocks_species.extend(block['seq_id'].map(lambda x: x[0:4]))

        # Compute features
        with mp.Pool(processes=num_processes) as pool:
            results = pd.DataFrame(pool.imap(seqfeat.feat_all, blocks_seq, chunksize=50))

        # Merge subsets and save
        results.set_index([blocks_id, blocks_class, blocks_species], inplace=True)
        results.index.names = ['block_id', 'ordered', 'species_id']
        results.to_csv(f'features_{i}.tsv', sep='\t')
