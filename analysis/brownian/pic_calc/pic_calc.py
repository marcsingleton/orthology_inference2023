"""Calculate phylogenetic contrasts from features."""

import Bio.Phylo as Phylo
import multiprocessing as mp
import pandas as pd
from copy import deepcopy
from os import environ


def get_contrasts(node, results):
    leaf1, leaf2 = node.clades

    if leaf1.is_terminal() and leaf2.is_terminal():
        val1, bl1 = leaf1.value, leaf1.branch_length
        val2, bl2 = leaf2.value, leaf2.branch_length
    elif leaf1.is_terminal():
        val1, bl1 = leaf1.value, leaf1.branch_length
        val2, bl2 = get_contrasts(leaf2, results)
    elif leaf2.is_terminal():
        val1, bl1 = get_contrasts(leaf1, results)
        val2, bl2 = leaf2.value, leaf2.branch_length
    else:
        val1, bl1 = get_contrasts(leaf1, results)
        val2, bl2 = get_contrasts(leaf2, results)

    bl_sum = bl1 + bl2
    value = (val1 * bl2 + val2 * bl1) / bl_sum
    branch_length = node.branch_length + bl1 * bl2 / bl_sum
    results.append((val1 - val2) / bl_sum)

    return value, branch_length


def block_contrasts(grouptree):
    # Initialize variables
    (key, block), tree = grouptree
    feature_pics = {}  # Dictionary of feature:PICs pairs
    terminals = {taxon_ids[terminal.name]: terminal for terminal in tree.get_terminals()}

    for feature in block:
        # Set terminal values on tree
        for species_id, value in zip(block.index.get_level_values('species_id'), block[feature]):
            terminals[species_id].value = value

        # Get PICs and store in dictionary
        pics = []
        get_contrasts(tree.root, pics)  # Pass root since tree object has no clades attribute
        feature_pics[feature] = pics

    # Convert dictionary to dataframe
    idx = [block.index.get_level_values('block_id')[:-1],
           block.index.get_level_values('ordered')[:-1],
           pd.Index([min(block.index.get_level_values('length'))] * (len(block) - 1), name='min_length')]
    return pd.DataFrame(feature_pics, index=idx)


# Input variables
path_features = '../feature_calc/features.tsv'
path_tree = '../prune_tree_25/out/drosophila-10spec-tree.nwk'
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

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    # Load data
    features = pd.read_csv(path_features, sep='\t', index_col=list(range(5)))
    tree = Phylo.read(path_tree, 'newick')

    # Compute PICs for each block
    groups = [group for group in features.groupby(level='block_id') if len(group[1]) == 10]  # Second item in tuple is df
    trees = [deepcopy(tree) for i in range(len(groups))]  # Deepcopy necessary due to nested structure of trees
    with mp.Pool(processes=num_processes) as pool:
        blocks_pics = pool.imap(block_contrasts, zip(groups, trees), chunksize=50)  # Dataframe of PICs for each block
        pics = pd.concat(blocks_pics)

    # Concatenate block dataframes and save
    pics.to_csv('pics.tsv', sep='\t')

"""
NOTES
The somewhat awkward zipping of the group and trees lists is due to a limitation of groupby and mp.Pool
    It is possible to filter groups with a built-in method, which is likely more efficient, but the filtration combines the data again, which would force a second round of grouping.
    Ideally, the code would pass the groups and the block_contrasts function would instantiate a new tree as needed.
        However, mp.Pool cannot parallelize functions which are not defined at initialization.
        Thus, I cannot create a block_contrasts factory function which can include the tree in the enclosing scope.
            The options are either to pass the tree as a argument or make block_contrasts dependent on the global environment.

DEPENDENCIES
../feature_calc/feature_calc.py
    ../feature_calc/features.tsv
../prune_tree_25/prune_tree_25.py
    ../prune_tree_25/out/drosophila-10spec-tree.nwk
"""