"""Calculate phylogenetic contrasts from raw features."""

import Bio.Phylo as Phylo
import multiprocessing as mp
import os
import pandas as pd
import re
from copy import deepcopy
from os import environ


def clade_contrasts(clade, pics):
    leaf1, leaf2 = clade.clades

    if leaf1.is_terminal() and leaf2.is_terminal():
        val1, bl1 = leaf1.value, leaf1.branch_length
        val2, bl2 = leaf2.value, leaf2.branch_length
    elif leaf1.is_terminal():
        val1, bl1 = leaf1.value, leaf1.branch_length
        val2, bl2 = clade_contrasts(leaf2, pics)
    elif leaf2.is_terminal():
        val1, bl1 = clade_contrasts(leaf1, pics)
        val2, bl2 = leaf2.value, leaf2.branch_length
    else:
        val1, bl1 = clade_contrasts(leaf1, pics)
        val2, bl2 = clade_contrasts(leaf2, pics)

    bl_sum = bl1 + bl2
    value = (val1 * bl2 + val2 * bl1) / bl_sum
    branch_length = clade.branch_length + bl1 * bl2 / bl_sum
    pics.append((val1 - val2) / bl_sum)

    return value, branch_length


def block_contrasts(grouproot):
    (key, block), root = grouproot
    feature_pics = {}  # Dictionary of feature:PICs pairs
    for feature in block:
        # Set tip values on tree
        for terminal in root.get_terminals():
            terminal.value = block.loc[(key[0], key[1], taxon_ids[terminal.name]), feature]

        # Get PICs and store in dictionary
        pics = []
        clade_contrasts(root, pics)
        feature_pics[feature] = pics

    # Convert dictionary to dataframe
    idx = block.index.droplevel('species_id')[:-1]  # Drop species_id and last row
    return pd.DataFrame(feature_pics, index=idx)


feature_dir = '../feature_calc_block/'  # Feature directory must end in /
num_processes = int(environ['SLURM_NTASKS'])
taxon_ids = {'ananassae': 7217,
             'erecta': 7220,
             'grimshawi': 7222,
             'melanogaster': 7227,
             'mojavensis': 7230,
             'persimilis': 7234,
             'pseudoobscura': 7237,
             'virilis': 7244,
             'yakuba': 7245,
             'willistoni': 7260}

tree = Phylo.read('../prune_d_tree_25/drosophila-10spec-tree.nwk', 'newick')
paths = filter(lambda x: re.match('features_[0-9]+\.tsv', x), os.listdir(feature_dir))
if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    for path in paths:
        # Load data
        df_raw = pd.read_csv(feature_dir + path, sep='\t', index_col=[0, 1, 2])
        df_raw.index.names = ['block_id', 'ordered', 'species_id']

        # Get file index
        j0 = path.find('_')
        j1 = path.find('.tsv')
        i = path[j0+1:j1]

        # Compute PICs for each block
        groups = df_raw.groupby(level=['block_id', 'ordered'])
        roots = [deepcopy(tree.clade) for i in range(len(groups))]  # Deepcopy necessary due to nested structure of trees
        with mp.Pool(processes=num_processes) as pool:
            blocks_pics = pool.imap(block_contrasts, zip(groups, roots), chunksize=50)  # Dataframe of PICs for each block
            df_pic = pd.concat(blocks_pics)

        # Concatenate block dataframes and save
        df_pic.to_csv(f'pics_{i}.tsv', sep='\t')
