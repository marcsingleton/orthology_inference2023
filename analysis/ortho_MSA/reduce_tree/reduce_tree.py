"""Reduce gene clusters to representative polypeptides using trees."""

import multiprocessing as mp
import os
import re
from itertools import combinations, product

import numpy as np
import skbio.stats.distance as distance
import skbio.tree
from src.utils import read_fasta


def get_ktuples(seq, k):
    """Return k-tuples of sequence."""
    ktuples = {}
    for i in range(len(seq)-k+1):
        ktuple = seq[i:i+k]
        ktuples[ktuple] = ktuples.get(ktuple, 0) + 1
    return ktuples


def get_ktuple_distance(seq1, seq2, k, p=1):
    """Return k-tuple distance between two sequences."""
    ktuples1 = get_ktuples(seq1, k)
    ktuples2 = get_ktuples(seq2, k)
    ktuples = set(ktuples1) | set(ktuples2)

    d = 0
    for ktuple in ktuples:
        d += abs(ktuples1.get(ktuple, 0) - ktuples2.get(ktuple, 0)) ** p
    return d


def get_max_gnid_distances(dm):
    """Return matrix of maximum distances between GNIDs."""
    gnid2idx = {}
    name2idx = {}
    for tip_name in dm.ids:
        gnid = tip_name.split(':')[1]
        try:
            idx = gnid2idx[gnid]
        except KeyError:
            idx = len(gnid2idx)
            gnid2idx[gnid] = idx
        name2idx[tip_name] = idx

    data = np.zeros((len(gnid2idx), len(gnid2idx)))
    for tip_name1, tip_name2 in product(dm.ids, repeat=2):
        idx1 = (name2idx[tip_name1], name2idx[tip_name2])
        idx2 = (tip_name1, tip_name2)
        data[idx1] = max(data[idx1], dm[idx2])
    np.fill_diagonal(data, 0)  # Set diagonal to 0
    return distance.DistanceMatrix(data, ids=sorted(gnid2idx, key=lambda x: gnid2idx[x]))  # Ensure sort


def update_max_gnid_distances(msd, dm, tree, gnids):
    """Update pairs containing elements of gnids in matrix of maximum distances."""
    data = -msd.data
    tip_names1 = [tip.name for tip in tree.tips()]
    tip_names2 = [tip_name for tip_name in tip_names1 if tip_name.split(':')[1] in gnids]
    name2idx = {tip_name: msd.index(tip_name.split(':')[1]) for tip_name in tip_names1}
    for tip_name1, tip_name2 in product(tip_names1, tip_names2):
        idx1 = (name2idx[tip_name1], name2idx[tip_name2])
        idx2 = (tip_name1, tip_name2)
        data[idx1] = max(data[idx1], dm[idx2])
        data[idx1[::-1]] = max(data[idx1], dm[idx2])
    np.fill_diagonal(data, 0)  # Set diagonal to 0
    data = np.abs(data)
    return distance.DistanceMatrix(data, ids=sorted(msd.ids, key=lambda x: msd.index(x)))  # Ensure sort


def update_tip_names(tree):
    """Cache tip names in internal nodes of tree."""
    for node in tree.postorder():
        if node.is_tip():
            node.tip_names = {node.name}
        else:
            node.tip_names = set().union(*[child.tip_names for child in node.children])


def reduce(OGid, OG):
    """Return representative PPIDs for all GNIDs associated with tips in node."""
    # Extract sequences
    seqs = []
    ids = []
    for ppid in OG:
        seqs.append(ppid2seq[ppid])
        gnid, spid = ppid2data[ppid]
        ids.append(f"'{spid}:{gnid}:{ppid}'")  # Wrap in quotes to ensure correct parsing

    # Make distance matrix
    k, p = 4, 2  # Tuple size and power
    i, j = 0, 1  # Matrix indices
    dm0 = np.zeros((len(seqs), len(seqs)))
    for seq1, seq2 in combinations(seqs, 2):
        # Calculate distance and store in matrix
        d = get_ktuple_distance(seq1.translate(table), seq2.translate(table), k, p)
        dm0[i, j] = d
        dm0[j, i] = d

        # Calculate indices
        j += 1
        if j > len(seqs) - 1:
            i += 1
            j = i + 1
    dm0 = distance.DistanceMatrix(dm0, ids=ids)

    # Make tree
    tree = skbio.tree.nj(dm0)
    update_tip_names(tree)
    dm = tree.tip_tip_distances()  # Use tree distances for pruning
    msd = get_max_gnid_distances(dm)

    # Prune tree
    gnids = {ppid2data[ppid][0] for ppid in OG}
    while len(tree.tip_names) > len(gnids):
        # Remove non-minimal tips in single-species clades
        for node in tree.postorder():
            if node.is_tip():
                node.min_tips = {(node.name, node.length)}
            elif len({name.split(':')[1] for child in node.children for name, _ in child.min_tips}) == 1:
                name, length = min([min_tip for child in node.children for min_tip in child.min_tips], key=lambda x: x[1])
                node.min_tips = {(name, length + node.length)}
            else:
                node.min_tips = {min_tip for child in node.children for min_tip in child.min_tips}
        min_names = {tip_name for tip_name, _ in tree.min_tips}
        if min_names < tree.tip_names:
            tip_gnids = {tip_name.split(':')[1] for tip_name in (tree.tip_names - min_names)}
            tree = tree.shear(min_names)
            update_tip_names(tree)
            msd = update_max_gnid_distances(msd, dm, tree, tip_gnids)

        # Split tree
        trees = []
        for node in tree.traverse(include_self=False):
            tip_gnids1 = {tip_name.split(':')[1] for tip_name in node.tip_names}
            tip_gnids2 = {tip_name.split(':')[1] for tip_name in (tree.tip_names - node.tip_names)}

            if tip_gnids1 == gnids:
                tree1 = tree.shear(node.tip_names)
                msd1 = update_max_gnid_distances(msd, dm, tree1, tip_gnids2)
                trees.append((tree1, msd1))
            if tip_gnids2 == gnids:
                tree2 = tree.shear(tree.tip_names - node.tip_names)
                msd2 = update_max_gnid_distances(msd, dm, tree2, tip_gnids1)
                trees.append((tree2, msd2))
        if trees:
            tree, msd = min(trees, key=lambda x: x[1].data.sum())
            update_tip_names(tree)

    # Extract sequences from tree
    reduced_OG = [tip_name.split(':') for tip_name in tree.tip_names]
    return OGid, reduced_OG


ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}
table = {ord('I'): '!', ord('L'): '!', ord('M'): '!', ord('V'): '!',
         ord('F'): '@', ord('W'): '@', ord('Y'): '@',
         ord('A'): '#', ord('S'): '#', ord('T'): '#',
         ord('D'): '$', ord('E'): '$', ord('N'): '$',
         ord('H'): '%', ord('K'): '%', ord('R'): '%', ord('Q'): '%',
         ord('C'): '^',
         ord('P'): '&',
         ord('G'): '~'}
num_processes = 2

# Load genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        genomes.append((spid, source, prot_path))

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, _ = line.split()
        ppid2data[ppid] = (gnid, spid)

# Load seqs
ppid2seq = {}
for _, source, prot_path in genomes:
    fasta = read_fasta(prot_path)
    for header, seq in fasta:
        ppid = re.search(ppid_regex[header], line).group(1)
        ppid2seq[ppid] = seq

# Load OGs
OGs = {}
with open('../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        _, OGid, _, edges = line.rstrip().split('\t')
        ppids = {node for edge in edges.split(',') for node in edge.split(':')}
        OGs[OGid] = ppids

if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        reduced_OGs = pool.starmap(reduce, OGs.items())

    # Write reduced clusters to file
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open('out/clusters.tsv', 'w') as file:
        file.write('OGid\tppids\n')
        for OGid, reduced_OG in reduced_OGs:
            nodestring = ','.join([record[2] for record in reduced_OG])
            file.write(f'{OGid}\t{nodestring}\n')

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/fasta/dmel-all-translation-r6.38.fasta
../../ortho_cluster3/cluster4+_graph/cluster4+_graph.py
    ../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../config/genomes.tsv
"""