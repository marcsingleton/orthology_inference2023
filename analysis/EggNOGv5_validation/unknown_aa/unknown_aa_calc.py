"""
Generate CSV of number of sequences with unknown amino acids and their IDs for each alignment.
Generate CSV of number, fraction, and contiguous lengths of unknown amino acids for each sequence.
"""

import gzip
import os
import pandas as pd
from Bio import AlignIO
from itertools import groupby

pathdir = [('../filter_count/out/7214_members/equal_+5_members.tsv', '../../../data/EggNOGv5/drosophilidae/7214/'),
           ('../../../data/EggNOGv5/drosophilidae/7214_members.tsv', '../../../data/EggNOGv5/drosophilidae/7214/')]

for path_members, dir_msa in pathdir:
    with open(path_members) as file_member:
        ali_data = []
        seq_data = []
        for line in file_member:
            fields = line.rstrip().split('\t')
            # Alignment metrics
            ali_id = fields[1]
            seqs = []  # Sequence IDs in alignment
            seqs_x = []  # Sequence IDs in alignment that contain unknown amino acids
            numseq = 0  # Number of sequences in alignment
            numseq_x = 0  # Number of sequences that contain unknown amino acids in alignment

            with gzip.open(dir_msa + ali_id + '.raw_alg.faa.gz', 'rt') as file_MSA:  # read, text mode
                MSA = AlignIO.read(file_MSA, 'fasta')
                for record in MSA:
                    # Sequence metrics
                    seq = str(record.seq.upper()).translate({ord('-'): None})  # Sequence with gaps removed and uppercase
                    seq_id = record.id
                    seq_len = len(seq)
                    num_x = 0  # Number of unknown amino acids in sequence
                    frac_x = 0  # Fraction of unknown amino acids in sequence
                    lens_x = []  # Lengths of contiguous segments of unknown amino acids

                    if 'X' in seq:
                        # Update X sequence metrics
                        num_x = seq.count('X')
                        frac_x = num_x / seq_len
                        groups_x = [list(g) for k, g in groupby(seq)]
                        for group_x in groups_x:
                            if set(group_x) == {'X'}:
                                lens_x.append(len(group_x))

                        # Update X alignment metrics
                        numseq_x += 1
                        seqs_x.append(seq_id)

                    # Update overall alignment metrics
                    numseq += 1
                    seqs.append(seq_id)

                    # Append current sequence metrics to list
                    seq_data.append({'seq_ID': seq_id, 'seq_len': seq_len, 'ali_ID': ali_id,
                                     'num_x': num_x, 'frac_x': frac_x, 'lens_x': lens_x})

            # Append current alignment metrics to list
            ali_data.append({'ali_id': ali_id, 'numseq': numseq, 'seqs': seqs, 'numseq_x': numseq_x, 'seqs_x': seqs_x})

    # Instantiate dataframes from lists
    ali_df = pd.DataFrame(ali_data)
    seq_df = pd.DataFrame(seq_data)

    # Save dataframes
    root, _ = os.path.splitext(os.path.basename(path_members))  # Get name of member file without extension
    if not os.path.exists('out/' + root):
        os.makedirs('out/' + root)  # Recursive folder creation

    ali_df.to_csv('out/' + root + '/ali_df.csv', index=False)
    seq_df.to_csv('out/' + root + '/seq_df.csv', index=False)
    ali_df.to_pickle('out/' + root + '/ali_df.pkl')
    seq_df.to_pickle('out/' + root + '/seq_df.pkl')

"""
DEPENDENCIES
../filter_count/filter_count.py
    ../filter_count/out/7214_members/equal_+5_members.tsv
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
    ../../../data/EggNOGv5/drosophilidae/7214/
"""