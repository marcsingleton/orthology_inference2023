"""Calculate contrasts using segments which pass quality filters."""

import multiprocessing as mp
import os

import pandas as pd
import skbio


def get_args(grouped, tree):
    for name, group in grouped:
        yield name, group, tree.copy()


def get_contrasts(node):
    child1, child2 = node.children
    if child1.is_tip():
        contrasts1, val1, bl1 = [], child1.value, child1.length
    else:
        contrasts1, val1, bl1 = get_contrasts(child1)
    if child2.is_tip():
        contrasts2, val2, bl2 = [], child2.value, child2.length
    else:
        contrasts2, val2, bl2 = get_contrasts(child2)

    bl_sum = bl1 + bl2
    value = (val1 * bl2 + val2 * bl1) / bl_sum
    branch_length = node.length + bl1 * bl2 / bl_sum
    contrasts = contrasts1 + contrasts2
    contrasts.append((val1 - val2) / bl_sum)

    return contrasts, value, branch_length


def apply_contrasts(args):
    name, group, tree = args

    # Map features to tips
    tree = tree.shear(group['spid'])
    spid2idx = {spid: i for i, spid in enumerate(group['spid'])}
    for tip in tree.tips():
        idx = spid2idx[tip.name]
        tip.value = group.iloc[idx].drop(['OGid', 'start', 'stop', 'ppid', 'spid'])

    # Get contrasts
    tree.length = 0  # Set root length to 0 for convenience
    contrasts, value, _ = get_contrasts(tree)

    return name, contrasts, value


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)

if __name__ == '__main__':
    # Load sequence data
    ppid2spid = {}
    with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            ppid, _, spid, _ = line.rstrip('\n').split('\t')
            ppid2spid[ppid] = spid

    # Load features
    features = pd.read_table('../get_features/out/features.tsv')
    features.loc[features['kappa'] == -1, 'kappa'] = 1
    features.loc[features['omega'] == -1, 'omega'] = 1
    feature_labels = list(features.columns.drop(['OGid', 'start', 'stop', 'ppid']))

    # Load segments
    rows = []
    with open('../aucpred_filter/out/regions_30.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            OGid, start, stop, disorder, ppids = line.rstrip('\n').split('\t')
            for ppid in ppids.split(','):
                rows.append({'OGid': OGid, 'start': int(start), 'stop': int(stop), 'ppid': ppid, 'spid': ppid2spid[ppid]})
    segments = pd.DataFrame(rows).merge(features, how='left', on=['OGid', 'start', 'stop', 'ppid'])
    regions = segments.groupby(['OGid', 'start', 'stop'])

    # Apply contrasts
    args = get_args(regions, tree_template)
    with mp.Pool(processes=num_processes) as pool:
        # Using imap distributes the construction of the args tuples
        # However to force computation before pool is closed we must call list on it
        records = list(pool.imap(apply_contrasts, args, chunksize=50))

    # Write contrasts and means to file
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open('out/contrasts.tsv', 'w') as file1, open('out/roots.tsv', 'w') as file2:
        file1.write('\t'.join(['OGid', 'start', 'stop', 'contrast_id'] + feature_labels) + '\n')
        file2.write('\t'.join(['OGid', 'start', 'stop'] + feature_labels) + '\n')
        for name, contrasts, value in records:
            # Write contrasts
            for idx, contrast in enumerate(contrasts):
                fields1 = [name[0], str(name[1]), str(name[2]), str(idx)]
                for feature_label in feature_labels:
                    fields1.append(str(contrast[feature_label]))
                file1.write('\t'.join(fields1) + '\n')

            # Write means
            fields2 = [name[0], str(name[1]), str(name[2])]
            for feature_label in feature_labels:
                fields2.append(str(value[feature_label]))
            file2.write('\t'.join(fields2) + '\n')

"""
DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../aucpred_filter/aucpred_filter.py
    ../aucpred_filter/out/regions_30.tsv
../get_features/get_features.py
    ../get_features/out/features.tsv
"""