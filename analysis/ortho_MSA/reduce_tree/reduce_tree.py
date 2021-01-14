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


def get_sum(tree, spids, dm):
    spid2idx = {spid: i for i, spid in enumerate(spids)}
    dm_spid = np.zeros((len(spid2idx), len(spid2idx)))
    for tip1, tip2 in product(tree.tips(), repeat=2):
        idx1 = (spid2idx[tip1.name.split(':')[0]], spid2idx[tip2.name.split(':')[0]])
        idx2 = (f"'{tip1.name}'", f"'{tip2.name}'")  # Wrap in quotes for consistency with dm
        dm_spid[idx1] = max(dm_spid[idx1], dm[idx2])
    np.fill_diagonal(dm_spid, 0)  # Set diagonal to 0
    return dm_spid.sum()


def reduce(OGid, OG):
    seqs = [seq for gnid in OG for seq in gnid2seqs[gnid]]
    spids = set([seq[1] for seq in seqs])

    # Make tree
    k, p = 4, 2  # Tuple size and power
    i, j = 0, 1  # Matrix indices
    dm = np.zeros((len(seqs), len(seqs)))
    for (ppid1, _, seq1), (ppid2, _, seq2) in combinations(seqs, 2):
        # Calculate distance and store in matrix
        d = get_ktuple_distance(seq1.translate(table), seq2.translate(table), k, p)
        dm[i, j] = d
        dm[j, i] = d

        # Calculate indices
        j += 1
        if j > len(seqs) - 1:
            i += 1
            j = i + 1
    dm = distance.DistanceMatrix(dm, ids=[f"'{spid}:{ppid}'" for ppid, spid, _ in seqs])
    tree = skbio.tree.nj(dm)

    # Split tree into subtrees
    while len(list(tree.tips())) > len(spids):
        # Mark node with tips
        for node in tree.traverse():
            node.tip_names = set([tip.name for tip in node.tips(include_self=True)])

        # Record subtrees
        subtrees = []
        for node in tree.traverse(include_self=False):
            tip_spids1 = set([tip_name.split(':')[0] for tip_name in node.tip_names])
            tip_spids2 = set([tip_name.split(':')[0] for tip_name in (tree.tip_names - node.tip_names)])

            if tip_spids2 == spids:
                subtree = tree.shear(tree.tip_names - node.tip_names)
                subtrees.append(subtree)
            if tip_spids1 == spids:
                subtree = tree.shear(node.tip_names)
                subtrees.append(subtree)
        tree = min(subtrees, key=lambda t: get_sum(t, spids, dm))

    # Extract sequences from tree
    rOG = []
    for tip in tree.tips():
        spid, ppid = tip.name.split(':')
        rOG.append((spid, ppid2gnid[ppid], ppid))
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

# Load pp metadata
ppid2gnid = {}
with open('../../ortho_search/ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        ppid, gnid, _ = line.split()
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
with open('../../ortho_cluster3/clique5+_community/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        OGs[OGid] = gnids

if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        rOGs = pool.starmap(reduce, OGs.items())

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
../../ortho_cluster3/clique5+_community/clique5+_community2.py
    ../../ortho_cluster3/clique5+_community/out/ggraph2/5clique/gclusters.txt
../../ortho_search/ppid2meta/ppid2meta.py
    ../../ortho_search/ppid2meta/out/ppid2meta.tsv
../config/genomes.tsv
"""