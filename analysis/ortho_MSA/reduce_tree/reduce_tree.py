"""Reduce gene clusters to representative polypeptides using trees."""

import multiprocessing as mp
import os
import re
from itertools import combinations, product

import numpy as np
import skbio.stats.distance as distance
import skbio.tree


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
            node.tip_names = set([node.name])
        else:
            node.tip_names = set().union(*[child.tip_names for child in node.children])


def reduce(OGid, OG):
    """Return representative PPIDs for all GNIDs associated with tips in node."""
    # Make distance matrix
    k, p = 4, 2  # Tuple size and power
    i, j = 0, 1  # Matrix indices
    seqs = [(ppid, gnid, spid, seq) for gnid in OG for ppid, spid, seq in gnid2seqs[gnid]]
    dm0 = np.zeros((len(seqs), len(seqs)))
    for (ppid1, _, _, seq1), (ppid2, _, _, seq2) in combinations(seqs, 2):
        # Calculate distance and store in matrix
        d = get_ktuple_distance(seq1.translate(table), seq2.translate(table), k, p)
        dm0[i, j] = d
        dm0[j, i] = d

        # Calculate indices
        j += 1
        if j > len(seqs) - 1:
            i += 1
            j = i + 1
    dm0 = distance.DistanceMatrix(dm0, ids=[f"'{spid}:{gnid}:{ppid}'" for ppid, gnid, spid, _ in seqs])  # Wrap in quotes to ensure correct parsing

    # Make tree
    tree = skbio.tree.nj(dm0)
    update_tip_names(tree)
    dm = tree.tip_tip_distances()  # Use tree distances for pruning
    msd = get_max_gnid_distances(dm)

    # Prune tree
    while len(tree.tip_names) > len(OG):
        # Remove non-minimal tips in single-species clades
        for node in tree.postorder():
            if node.is_tip():
                node.min_tips = set([(node.name, node.length)])
            elif len(set([name.split(':')[1] for child in node.children for name, _ in child.min_tips])) == 1:
                name, length = min([min_tip for child in node.children for min_tip in child.min_tips], key=lambda x: x[1])
                node.min_tips = set([(name, length + node.length)])
            else:
                node.min_tips = set([min_tip for child in node.children for min_tip in child.min_tips])
        min_names = set([tip_name for tip_name, _ in tree.min_tips])
        if min_names < tree.tip_names:
            tip_gnids = set([tip_name.split(':')[1] for tip_name in (tree.tip_names - min_names)])
            tree = tree.shear(min_names)
            update_tip_names(tree)
            msd = update_max_gnid_distances(msd, dm, tree, tip_gnids)

        # Split tree
        trees = []
        for node in tree.traverse(include_self=False):
            tip_gnids1 = set([tip_name.split(':')[1] for tip_name in node.tip_names])
            tip_gnids2 = set([tip_name.split(':')[1] for tip_name in (tree.tip_names - node.tip_names)])

            if tip_gnids1 == OG:
                tree1 = tree.shear(node.tip_names)
                msd1 = update_max_gnid_distances(msd, dm, tree1, tip_gnids2)
                trees.append((tree1, msd1))
            if tip_gnids2 == OG:
                tree2 = tree.shear(tree.tip_names - node.tip_names)
                msd2 = update_max_gnid_distances(msd, dm, tree2, tip_gnids1)
                trees.append((tree2, msd2))
        if trees:
            tree, msd = min(trees, key=lambda x: x[1].data.sum())
            update_tip_names(tree)

    # Extract sequences from tree
    rOG = [tip_name.split(':') for tip_name in tree.tip_names]
    return OGid, rOG


pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}
table = {ord('I'): '!', ord('L'): '!', ord('M'): '!', ord('V'): '!',
         ord('F'): '@', ord('W'): '@', ord('Y'): '@',
         ord('A'): '#', ord('S'): '#', ord('T'): '#',
         ord('D'): '$', ord('E'): '$', ord('N'): '$',
         ord('H'): '%', ord('K'): '%', ord('R'): '%', ord('Q'): '%',
         ord('C'): '^',
         ord('P'): '&',
         ord('G'): '~'}
num_processes = int(os.environ['SLURM_NTASKS'])

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        genomes.append((spid, source, prot_path))

# Load seq metadata
ppid2gnid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, _, _ = line.split()
        ppid2gnid[ppid] = gnid

# Load seqs
gnid2seqs = {}
for spid0, source, prot_path in genomes:
    with open(prot_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid0 = re.search(pp_regex[source], line).group(1)
                gnid = ppid2gnid[ppid0]
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq0 = ''.join(seqlines)

            try:
                for _, _, seq1 in gnid2seqs[gnid]:
                    if seq0 == seq1:
                        break
                else:
                    gnid2seqs[gnid].append((ppid0, spid0, seq0))
            except KeyError:
                gnid2seqs[gnid] = [(ppid0, spid0, seq0)]

# Load OGs
OGs = {}
with open('../../ortho_cluster3/clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        OGs[OGid] = gnids

if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        rOGs = pool.starmap(reduce, list(OGs.items()))

    # Make output directory
    if not os.path.exists('out/'):
        os.mkdir('out/')

    # Write reduced clusters to file
    with open('out/rclusters.tsv', 'w') as outfile:
        outfile.write('OGid\tspid\tgnid\tppid\n')
        for OGid, rOG in rOGs:
            for entry in rOG:
                outfile.write(OGid + '\t' + '\t'.join(entry) + '\n')

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_cluster3/clique4+_gcommunity/clique4+_gcommunity2.py
    ../../ortho_cluster3/clique4+_gcommunity/out/ggraph2/5clique/gclusters.txt
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../config/genomes.tsv
"""