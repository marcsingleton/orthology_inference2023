"""Simulate sequence evolution."""

import os
import re
from copy import deepcopy
from io import StringIO

import numpy as np
import scipy.stats as stats
import skbio
from scipy.special import gammainc
from scipy.stats import gamma
from src.evosim.asr import get_tree
from src.utils import read_fasta


class SeqEvolver:
    """A class for mutating sequences.

    Insertions and deletions are right-oriented. For insertions, the residues
    are inserted immediately left of the given index. For deletions, active
    residues are deleted beginning at the given index and moving left until the
    number of residues deleted is equal to the given length or no more residues
    remain.

    Parameters
    ----------
    seq: ndarray
        One-dimensional array where symbols are stored as numerical indices.
        The indices must match the order given in jump_matrices and sym_dists.
    rates: ndarray
        Array with shape (3, len(seq)). The rows correspond to substitution,
        insertion, and deletion rates, respectively. The columns correspond to
        the index in seq.
    activities: ndarray
        One-dimenionsal array with boolean values indicating if residue is
        active. Deleted residues are inactive.
    residue_ids: ndarray
        One-dimensional array with unique integer identifier for each symbol in
        seq.
    partition_ids: ndarray
        One-dimenionsal array with the partition integer identifier for each
        symbol in seq.
    jump_matrices: dict of ndarrays
        Dict keyed by partition_id where values are jump distributions.
    sym_dists: dict of nd arrays
        Dict keyed by partition_id where values are symbol distributions.
    insertion_dists: dict of rv_discrete
        Dict keyed by partition_id where values are rv_discrete for generating
        random insertion lengths. Values must support a rvs method.
    deletion_dists: dict of rv_discrete
        Dict keyed by partition_id where values are rv_discrete for generating
        random deletion lengths. Values must support a rvs method.
    """
    def __init__(self, seq, rates, activities, residue_ids, partition_ids, jump_matrices, sym_dists, insertion_dists, deletion_dists):
        self.seq = seq
        self.rates = rates
        self.activities = activities
        self.residue_ids = residue_ids
        self.partition_ids = partition_ids
        self.jump_matrices = jump_matrices
        self.sym_dists = sym_dists
        self.insertion_dists = insertion_dists
        self.deletion_dists = deletion_dists

    def __deepcopy__(self, memodict={}):
        seq = np.copy(self.seq)
        rates = np.copy(self.rates)
        activities = np.copy(self.activities)
        residue_ids = np.copy(self.residue_ids)
        partition_ids = np.copy(self.partition_ids)
        return SeqEvolver(seq, rates, activities, residue_ids, partition_ids, self.jump_matrices, self.sym_dists, self.insertion_dists, self.deletion_dists)

    def mutate(self, residue_index):
        """Mutate the sequence."""
        event_ids = np.arange(3*len(self.seq))
        active_rates = (self.rates * self.activities).flatten()
        p = active_rates / active_rates.sum()
        event_id = rng.choice(event_ids, p=p)
        i, j = event_id // len(self.seq), event_id % len(self.seq)
        if i == 0:
            return self.substitute(j, residue_index)
        elif i == 1:
            return self.insert(j, residue_index)
        elif i == 2:
            return self.delete(j, residue_index)

    def substitute(self, j, residue_index):
        """Substitute residue at index j."""
        jump_dist = self.jump_matrices[self.partition_ids[j]][self.seq[j]]
        self.seq[j] = rng.choice(np.arange(len(jump_dist)), p=jump_dist)

        return residue_index

    def insert(self, j, residue_index):
        """Insert randomly generated residues at an index j."""
        partition_id = self.partition_ids[j]
        sym_dist = self.sym_dists[partition_id]
        length = self.insertion_dists[partition_id].rvs()

        # Make insertion arrays
        seq = rng.choice(np.arange(len(sym_dist)), size=length, p=sym_dist)
        rates = np.full((3, length), np.expand_dims(self.rates[:, j], axis=-1))
        activities = np.full(length, True)
        residue_ids = np.arange(residue_index, residue_index+length)
        partition_ids = np.full(length, partition_id)

        # Insert insertion into arrays
        self.seq = np.insert(self.seq, j+1, seq)
        self.rates = np.insert(self.rates, [j+1], rates, axis=1)  # Array insertion requires some special syntax
        self.activities = np.insert(self.activities, j+1, activities)
        self.residue_ids = np.insert(self.residue_ids, j+1, residue_ids)
        self.partition_ids = np.insert(self.partition_ids, j+1, partition_ids)

        return residue_index + length

    def delete(self, j, residue_index):
        """Delete residues at an index j."""
        partition_id = self.partition_ids[j]
        length = self.deletion_dists[partition_id].rvs()

        while length > 0 and j < len(self.seq):
            if self.activities[j]:
                self.activities[j] = False
                length -= 1
            j += 1

        return residue_index


def simulate_branch(node, evoseq, residue_index, t=None):
    if t is None:
        scale = 1/(evoseq.rates*evoseq.activities).sum()
        t = stats.expon.rvs(scale=scale)
    while t < node.length:
        residue_index = evoseq.mutate(residue_index)
        scale = 1/(evoseq.rates*evoseq.activities).sum()
        t += stats.expon.rvs(scale=scale)
    if node.children:
        if rng.random() > 0.5:
            child1, child2 = node.children
        else:
            child2, child1 = node.children
        evoseq1, evoseq2 = evoseq, deepcopy(evoseq)  # One sequence is the original, the other is copied
        residue_index, evoseqs1 = simulate_branch(child1, evoseq1, residue_index, t=t-node.length)
        residue_index, evoseqs2 = simulate_branch(child2, evoseq2, residue_index)
        evoseqs = evoseqs1 + evoseqs2
    else:
        evoseqs = [(node.name, evoseq)]

    return residue_index, evoseqs


def load_model(path):
    with open(path) as file:
        # Parse exchangeability matrix
        matrix = np.zeros((len(alphabet), len(alphabet)))
        for i in range(19):
            line = file.readline()
            for j, value in enumerate(line.split()):
                matrix[i + 1, j] = float(value)
                matrix[j, i + 1] = float(value)

        # Parse equilibrium frequencies
        for _ in range(2):
            line = file.readline()
        freqs = np.array([float(value) for value in line.split()])
        freqs = freqs / freqs.sum()  # Re-normalize to ensure sums to 1 (in case of floating point errors)
    matrix = freqs * matrix
    matrix = matrix / np.expand_dims(matrix.sum(axis=1), -1)  # Normalize rows to obtain jump matrix
    return matrix, freqs


rng = np.random.default_rng()
alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
sym2idx = {sym: i for i, sym in enumerate(alphabet)}
idx2sym = {i: sym for i, sym in enumerate(alphabet)}
insertion_rate = 0.008
deletion_rate = 0.01

# Load models
WAG_matrix, WAG_freqs = load_model('../config/WAG.txt')
disorder_matrix, disorder_freqs = load_model('../config/50red_D.txt')
jump_matrices = {1: WAG_matrix, 2: disorder_matrix}
sym_dists = {1: WAG_freqs, 2: disorder_freqs}

# Make indel dists
insertion_dists = {1: stats.geom(0.8), 2: stats.geom(0.75)}
deletion_dists = {1: stats.geom(0.6), 2: stats.geom(0.55)}

if not os.path.exists('out/'):
    os.mkdir('out/')

for path in os.listdir('../asr_generate/out/'):
    OGid = path[:4]
    fasta = read_fasta(f'../asr_generate/out/{OGid}_sample.mfa')
    length = len(fasta[0][1])

    # Load trees
    tree1 = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
    with open(f'../asr_indel/out/{OGid}.iqtree') as file:
        line = file.readline()
        while line != 'Tree in newick format:\n':
            line = file.readline()
        for _ in range(2):
            line = file.readline()
    tree2 = skbio.read(StringIO(line), 'newick', skbio.TreeNode)
    tree = get_tree(tree1, tree2)
    tree.length = 0  # Set root branch length to 0

    # Load partition model parameters
    partitions = {}
    with open(f'../asr_aa/out/{OGid}.iqtree') as file:
        # Get partition ID and name
        line = file.readline()
        while line.split() != ['ID', 'Name', 'Type', 'Seq', 'Site', 'Unique', 'Infor', 'Invar', 'Const']:  # Spacing can differ so check for fields
            line = file.readline()
        line = file.readline()
        while line != '\n':
            fields = line.split()  # File is not explicitly delimited, so just split on whitespace
            partition_id, name = int(fields[0]), fields[1]
            partitions[partition_id] = {'name': name}
            line = file.readline()

        # Get partition model parameters
        while line != '  ID  Model           Speed  Parameters\n':
            line = file.readline()
        line = file.readline()
        while line != '\n':
            fields = line.split()  # File is not explicitly delimited, so just split on whitespace
            partition_id, speed, parameters = int(fields[0]), float(fields[2]), fields[3]
            match = re.search(r'(?P<model>[^+]+)\+I{(?P<pinv>[0-9.e-]+)}\+G(?P<num_categories>[0-9]+){(?P<alpha>[0-9.e-]+)}', parameters)
            partition = partitions[partition_id]
            partition.update({'model': match['model'], 'speed': speed,
                              'pinv': float(match['pinv']), 'alpha': float(match['alpha']), 'num_categories': int(match['num_categories'])})
            line = file.readline()

    # Load partition regions
    partition_ids = np.empty(length)
    with open(f'../asr_aa/out/{OGid}.nex') as file:
        partition_id = 1
        for line in file:
            if 'charset' in line:
                match = re.search(r'charset (?P<name>[a-zA-Z0-9]+) = (?P<regions>[0-9 -]+);', line)
                for region in match['regions'].split():
                    start, stop = region.split('-')
                    partition_ids[int(start)-1:int(stop)] = partition_id
                partition_id += 1

    # Load rate categories
    for partition in partitions.values():
        pinv, alpha, num_categories = partition['pinv'], partition['alpha'], partition['num_categories']
        igfs = []  # Incomplete gamma function evaluations
        for i in range(num_categories+1):
            x = gamma.ppf(i/num_categories, a=alpha, scale=1/alpha)
            igfs.append(gammainc(alpha+1, alpha*x))
        rates = [(0, pinv)]
        for i in range(num_categories):
            rate = num_categories/(1-pinv) * (igfs[i+1] - igfs[i])
            rates.append((partition['speed'] * rate, (1-pinv)/num_categories))
        partition['rates'] = rates

    # Evolve sequences along tree
    residue_ids = np.arange(length)
    aa_dist = np.load(f'../asr_root/out/{OGid}_aa.npy')
    for header, seq in fasta:
        # Construct sequence object
        seq = np.array([sym2idx.get(sym, -1) for sym in seq])  # Use -1 for gap characters
        activities = np.array([False if sym == -1 else True for sym in seq])
        rates = np.empty((3, length))
        for j, (i, partition_id) in enumerate(zip(seq, partition_ids)):
            ps = aa_dist[:, i, j] / aa_dist[:, i, j].sum()  # Posterior for rate categories given symbol
            rs = np.array([r for r, _ in partitions[partition_id]['rates']])  # Rates of rate categories
            rate = (ps*rs).sum()
            rates[:, j] = [rate, insertion_rate*rate, deletion_rate*rate]
        evoseq = SeqEvolver(seq, rates, activities, residue_ids, partition_ids, jump_matrices, sym_dists, insertion_dists, deletion_dists)

        # Evolve! (and extract results)
        _, evoseqs = simulate_branch(tree, evoseq, length)

        unaligned_records = []
        for spid, evoseq in evoseqs:
            seq = []
            for idx, activity in zip(evoseq.seq, evoseq.activities):
                if activity:
                    seq.append(idx2sym.get(idx, '-'))
                else:
                    seq.append('-')
            unaligned_records.append((spid, seq, list(evoseq.residue_ids)))

        # Align sequences
        aligned_records = [(spid, []) for spid, _, _ in unaligned_records]
        spid2idx = {spid: 0 for spid, _, _ in unaligned_records}
        while any([spid2idx[spid] < len(seq) for spid, seq, _ in unaligned_records]):
            # Collect symbols and residue ids
            idx_records = []
            for spid, seq, residue_ids in unaligned_records:
                idx = spid2idx[spid]
                idx_records.append((seq[idx], residue_ids[idx], spid))

            # Find spids with priority id
            max_id = max([residue_id for _, residue_id, _ in idx_records])
            spids = {spid for _, residue_id, spid in idx_records if residue_id == max_id}

            # Append symbols to sequences
            for (spid, seq), (sym, _, _) in zip(aligned_records, idx_records):
                if spid in spids:
                    seq.append(sym)
                    spid2idx[spid] += 1
                else:
                    seq.append('-')

        # Write alignment
        with open(f'out/{OGid}_{header[1:]}.mfa', 'w') as file:
            for spid, seq in sorted(aligned_records):
                seqstring = '\n'.join([''.join(seq)[i:i+80] for i in range(0, len(seq), 80)]) + '\n'
                file.write(f'>{spid}\n' + seqstring)

"""
DEPENDENCIES
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../asr_generate/asr_generate.py
    ../asr_generate/out/*_sample.mfa
../asr_aa/asr_aa.py
    ../asr_aa/out/*.iqtree
    ../asr_aa/out/*.nex
../asr_root/aa.py
    ../asr_root/out/*_aa.npy
../config/50red_D.txt
../config/WAG.txt
"""