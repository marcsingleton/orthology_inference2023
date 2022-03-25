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

    An "immortal link" (Thorne et al., 1991) can be included at the beginning
    to allow indels before the first symbol. The immortal link is never
    inactive to ensure a sequence can always generate new symbols. Immortal
    links are not officially supported in this implementation, but they can be
    easily simulated by including an additional symbol at the beginning with a
    non-zero insertion rate and a deletion rate of 0. As a result of how the
    rate array is initialized, the symbol of the immortal link must have an
    associated rate in the rate matrix. However, to keep this rate constant,
    it is recommended to set the substitution rate to 0 as well.

    Parameters
    ----------
    seq: ndarray
        One-dimensional array where symbols are stored as numerical indices.
        The indices must match the order given in jump_matrices and sym_dists.
    rate_coefficients: ndarray
        Array with shape (3, len(seq)). The rows correspond to substitution,
        insertion, and deletion rate coefficients, respectively. The columns
        correspond to the index in seq.
    activities: ndarray
        One-dimensional array with boolean values indicating if the symbol is
        active. Deleted symbols are inactive.
    residue_ids: ndarray
        One-dimensional array with unique integer identifier for each symbol in
        seq.
    partition_ids: ndarray
        One-dimensional array with the partition integer identifier for each
        symbol in seq.
    rate_matrices: dict of ndarrays
        Dict keyed by partition_id where values are normalized rate matrices.
    sym_dists: dict of nd arrays
        Dict keyed by partition_id where values are symbol distributions.
    insertion_dists: dict of rv_discrete
        Dict keyed by partition_id where values are rv_discrete for generating
        random insertion lengths. Values must support an rvs method.
    deletion_dists: dict of rv_discrete
        Dict keyed by partition_id where values are rv_discrete for generating
        random deletion lengths. Values must support an rvs method.
    """
    def __init__(self, seq, rate_coefficients, activities, residue_ids, partition_ids,
                 rate_matrices, sym_dists, insertion_dists, deletion_dists):
        self.seq = seq
        self.rate_coefficients = rate_coefficients
        self.activities = activities
        self.residue_ids = residue_ids
        self.partition_ids = partition_ids
        self.rate_matrices = rate_matrices
        self.sym_dists = sym_dists
        self.insertion_dists = insertion_dists
        self.deletion_dists = deletion_dists

        # Calculate rates from rate matrices and rate coefficients
        rates = np.empty(len(seq))
        for j, (idx, partition_id) in enumerate(zip(self.seq, self.partition_ids)):
            rate = 0 if idx == -1 else -self.rate_matrices[partition_id][idx, idx]  # Check for "out-of-alphabet" symbols
            rates[j] = rate
        self.rates = rates * self.rate_coefficients

        # Calculate jump matrices from rate matrices
        jump_matrices = {}
        for partition_id, rate_matrix in rate_matrices.items():
            jump_matrix = np.copy(rate_matrix)
            np.fill_diagonal(jump_matrix, 0)
            jump_matrix = jump_matrix / np.expand_dims(jump_matrix.sum(axis=1), -1)  # Normalize rows to obtain jump matrix
            jump_matrices[partition_id] = jump_matrix
        self.jump_matrices = jump_matrices

    def __deepcopy__(self, memodict={}):
        seq = np.copy(self.seq)
        rate_coefficients = np.copy(self.rate_coefficients)
        activities = np.copy(self.activities)
        residue_ids = np.copy(self.residue_ids)
        partition_ids = np.copy(self.partition_ids)
        return SeqEvolver(seq, rate_coefficients, activities, residue_ids, partition_ids, self.rate_matrices, self.sym_dists, self.insertion_dists, self.deletion_dists)

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
        partition_id = self.partition_ids[j]
        jump_dist = self.jump_matrices[partition_id][self.seq[j]]
        idx = rng.choice(np.arange(len(jump_dist)), p=jump_dist)
        rate = -self.rate_matrices[partition_id][idx, idx]

        self.seq[j] = idx
        self.rates[:, j] = rate * self.rate_coefficients[:, j]

        return residue_index

    def insert(self, j, residue_index):
        """Insert randomly generated residues at an index j."""
        partition_id = self.partition_ids[j]
        sym_dist = self.sym_dists[partition_id]
        length = self.insertion_dists[partition_id].rvs()

        # Make insertion arrays
        # (rate_coefficients and rates are transposed since it makes the insertion syntax simpler)
        seq = rng.choice(np.arange(len(sym_dist)), size=length, p=sym_dist)
        rate_coefficients = np.full((length, 3), self.rate_coefficients[:, j])
        activities = np.full(length, True)
        residue_ids = np.arange(residue_index, residue_index+length)
        partition_ids = np.full(length, partition_id)
        rates = np.expand_dims([-self.rate_matrices[partition_id][idx, idx] for idx in seq], -1)

        # Insert insertion into arrays
        self.seq = np.insert(self.seq, j+1, seq)
        self.rate_coefficients = np.insert(self.rate_coefficients, j+1, rate_coefficients, axis=1)  # Insert matches first axis of inserted array to one given by axis keyword
        self.activities = np.insert(self.activities, j+1, activities)
        self.residue_ids = np.insert(self.residue_ids, j+1, residue_ids)
        self.partition_ids = np.insert(self.partition_ids, j+1, partition_ids)
        self.rates = np.insert(self.rates, j+1, rates * rate_coefficients, axis=1)

        return residue_index + length

    def delete(self, j, residue_index):
        """Delete residues at an index j."""
        partition_id = self.partition_ids[j]
        length = self.deletion_dists[partition_id].rvs()

        self.activities[j:j+length] = False

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
        for i in range(len(alphabet)-1):
            line = file.readline()
            for j, value in enumerate(line.split()):
                matrix[i + 1, j] = float(value)
                matrix[j, i + 1] = float(value)

        # Parse equilibrium frequencies
        for _ in range(2):
            line = file.readline()
        freqs = np.array([float(value) for value in line.split()])
    matrix = freqs * matrix
    rate = (freqs * matrix.sum(axis=1)).sum()
    matrix = matrix / rate  # Normalize average rate to 1
    np.fill_diagonal(matrix, -matrix.sum(axis=1))
    return matrix, freqs


rng = np.random.default_rng()
alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
sym2idx = {sym: i for i, sym in enumerate(alphabet)}
idx2sym = {i: sym for i, sym in enumerate(alphabet)}
insertion_rate = 0.008  # Indel rates are relative to the substitution rate scaling factor
deletion_rate = 0.01

# Load models
WAG_matrix, WAG_freqs = load_model('../config/WAG.txt')
disorder_matrix, disorder_freqs = load_model('../config/50red_D.txt')
rate_matrices = {1: WAG_matrix,
                 2: disorder_matrix}
sym_dists = {1: WAG_freqs / WAG_freqs.sum(),  # Re-normalize to ensure sums to 1
             2: disorder_freqs / disorder_freqs.sum()}

# Make indel dists
insertion_dists = {1: stats.geom(0.8), 2: stats.geom(0.75)}
deletion_dists = {1: stats.geom(0.6), 2: stats.geom(0.55)}

if not os.path.exists('out/'):
    os.mkdir('out/')

for path in [path for path in os.listdir('../asr_generate/out/') if path.endswith('_sample.mfa')]:
    # Load data and calculate "global" variables
    OGid = path.removesuffix('_sample.mfa')
    fasta = read_fasta(f'../asr_generate/out/{OGid}_sample.mfa')
    aa_dist = np.load(f'../asr_root/out/{OGid}_aa.npy')
    length = len(fasta[0][1])
    residue_ids = np.arange(-1, length)

    # Load trees
    tree1 = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
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
    partition_template = np.empty(length)
    with open(f'../asr_aa/out/{OGid}.nex') as file:
        partition_id = 1
        for line in file:
            if 'charset' in line:
                match = re.search(r'charset (?P<name>[a-zA-Z0-9]+) = (?P<regions>[0-9 -]+);', line)
                for region in match['regions'].split():
                    start, stop = region.split('-')
                    partition_template[int(start)-1:int(stop)] = partition_id
                partition_id += 1

    # Load rate categories
    # In IQ-TREE, only the shape parameter is fit and the rate parameter beta is set to alpha so the mean of gamma distribution is 1
    # The calculations here directly correspond to equation 10 in Yang. J Mol Evol (1994) 39:306-314.
    # Note the equation has a small typo where the difference in gamma function evaluations should be divided by the probability
    # of that category since technically it is the rate given that category
    # The order of arguments is also exchanged between the reference and scipy's implementation
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
    for header, seq in fasta:
        # Construct root-specific SeqEvolver arrays
        seq = np.array([sym2idx.get(sym, -1) for sym in seq])  # Use -1 for gap symbols
        activities = np.array([False if sym == -1 else True for sym in seq])
        rate_coefficients = np.empty((3, length))
        for j, (i, partition_id) in enumerate(zip(seq, partition_template)):
            ps = aa_dist[:, i, j] / aa_dist[:, i, j].sum()  # Posterior for rate categories given symbol
            rs = np.array([r for r, _ in partitions[partition_id]['rates']])  # Rates of rate categories
            rate = (ps*rs).sum()
            rate_coefficients[:, j] = [rate, insertion_rate*rate, deletion_rate*rate]

        # Insert immortal link
        j = np.nonzero(activities)[0][0]  # Index of first active symbol
        seq = np.insert(seq, 0, seq[j])
        rate_coefficients = np.insert(rate_coefficients, 0, [0, rate_coefficients[1, j], 0], axis=1)
        activities = np.insert(activities, 0, True)
        partition_ids = np.insert(partition_template, 0, partition_template[j])

        evoseq = SeqEvolver(seq, rate_coefficients, activities, residue_ids, partition_ids,
                            rate_matrices, sym_dists, insertion_dists, deletion_dists)

        # Evolve! (and extract results)
        _, evoseqs = simulate_branch(tree, evoseq, length)

        unaligned_records = []
        for spid, evoseq in evoseqs:
            seq = []
            for idx, activity in zip(evoseq.seq[1:], evoseq.activities[1:]):  # Exclude immortal link
                if activity:
                    seq.append(idx2sym.get(idx, '-'))
                else:
                    seq.append('-')
            unaligned_records.append((spid, seq, list(evoseq.residue_ids[1:])))

        # Align sequences
        aligned_ids = []
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
            aligned_ids.append(str(max_id))

            # Append symbols to sequences
            for (spid, seq), (sym, residue_id, _) in zip(aligned_records, idx_records):
                if spid in spids:
                    seq.append(sym)
                    spid2idx[spid] += 1
                else:
                    seq.append('-')

        # Write alignments
        with open(f'out/{OGid}_{header[1:]}.txt', 'w') as file:
            file.write(','.join(aligned_ids) + '\n')
        with open(f'out/{OGid}_{header[1:]}.mfa', 'w') as file:
            for spid, seq in sorted(aligned_records):
                seqstring = '\n'.join([''.join(seq)[i:i+80] for i in range(0, len(seq), 80)]) + '\n'
                file.write(f'>{spid}\n' + seqstring)

"""
DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
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