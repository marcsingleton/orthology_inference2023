"""Fit HMM to segmented alignments using conditional maximum likelihood via gradient descent."""

import json
import multiprocessing as mp
import os
import re
from functools import reduce

import numpy as np
import skbio
import src.brownian.hmm as hmm
import utils
from numpy import exp, log
from src.utils import read_fasta


# Gradient functions
def get_tree_prime_pi(tree, q0, q1):
    """Return derivative of probability of tree relative to pi."""
    s, conditional = utils.get_conditional(tree, q0, q1)
    d = ((exp(s) * conditional) * [[-1], [1]]).sum(axis=0)  # Broadcasting magic
    return d


def get_tree_prime_q0(tree, pi, q0, q1):
    """Return derivative of probability of tree relative to q0."""
    derivative = get_conditional_prime_q0(tree, q0, q1)
    d = (derivative * [[1-pi], [pi]]).sum(axis=0)  # Broadcasting magic
    return d


def get_tree_prime_q1(tree, pi, q0, q1):
    """Return derivative of probability of tree relative to q1."""
    derivative = get_conditional_prime_q1(tree, q0, q1)
    d = (derivative * [[1-pi], [pi]]).sum(axis=0)  # Broadcasting magic
    return d


def get_conditional_prime_q0(node, q0, q1):
    """Return derivative of conditional probabilities relative to q0."""
    child1, child2 = node.children

    if child1.is_tip():
        s1, conditional1 = 0, child1.conditional
        derivative1 = np.zeros((2, conditional1.shape[1]))
    else:
        s1, conditional1 = utils.get_conditional(child1, q0, q1)
        conditional1 = exp(s1) * conditional1  # Un-normalize
        derivative1 = get_conditional_prime_q0(child1, q0, q1)
    if child2.is_tip():
        s2, conditional2 = 0, child2.conditional
        derivative2 = np.zeros((2, conditional2.shape[1]))
    else:
        s2, conditional2 = utils.get_conditional(child2, q0, q1)
        conditional2 = exp(s2) * conditional2  # Un-normalize
        derivative2 = get_conditional_prime_q0(child2, q0, q1)

    p1 = utils.get_transition_matrix(q0, q1, child1.length)
    p2 = utils.get_transition_matrix(q0, q1, child2.length)
    d1 = get_transition_matrix_prime_q0(q0, q1, child1.length)
    d2 = get_transition_matrix_prime_q0(q0, q1, child2.length)
    term1 = np.matmul(d1, conditional1) + np.matmul(p1, derivative1)
    term2 = np.matmul(d2, conditional2) + np.matmul(p2, derivative2)
    derivative = term1*np.matmul(p2, conditional2) + np.matmul(p1, conditional1)*term2

    return derivative


def get_conditional_prime_q1(node, q0, q1):
    """Return derivative of conditional probabilities relative to q1."""
    child1, child2 = node.children

    if child1.is_tip():
        s1, conditional1 = 0, child1.conditional
        derivative1 = np.zeros((2, conditional1.shape[1]))
    else:
        s1, conditional1 = utils.get_conditional(child1, q0, q1)
        conditional1 = exp(s1) * conditional1  # Un-normalize
        derivative1 = get_conditional_prime_q1(child1, q0, q1)
    if child2.is_tip():
        s2, conditional2 = 0, child2.conditional
        derivative2 = np.zeros((2, conditional2.shape[1]))
    else:
        s2, conditional2 = utils.get_conditional(child2, q0, q1)
        conditional2 = exp(s2) * conditional2  # Un-normalize
        derivative2 = get_conditional_prime_q1(child2, q0, q1)

    p1 = utils.get_transition_matrix(q0, q1, child1.length)
    p2 = utils.get_transition_matrix(q0, q1, child2.length)
    d1 = get_transition_matrix_prime_q1(q0, q1, child1.length)
    d2 = get_transition_matrix_prime_q1(q0, q1, child2.length)
    term1 = np.matmul(d1, conditional1) + np.matmul(p1, derivative1)
    term2 = np.matmul(d2, conditional2) + np.matmul(p2, derivative2)
    derivative = term1*np.matmul(p2, conditional2) + np.matmul(p1, conditional1)*term2

    return derivative


def get_transition_matrix_prime_q0(q0, q1, t):
    """Return derivative transition matrix for two-state CTMC relative to q0."""
    q = q0 + q1
    d00 = -(q1 + (t*q0**2+t*q0*q1-q1)*exp(-q*t)) / q**2
    d01 = -d00
    d11 = q1*(1 - (q0*t+q1*t+1)*exp(-q*t)) / q**2
    d10 = -d11
    return np.array([[d00, d01], [d10, d11]])


def get_transition_matrix_prime_q1(q0, q1, t):
    """Return derivative transition matrix for two-state CTMC relative to q1."""
    q = q0 + q1
    d00 = q0*(1 - (t*q0+t*q1+1)*exp(-q*t)) / q**2
    d01 = -d00
    d11 = -(q0 + (t*q1**2+t*q0*q1-q0)*exp(-q*t)) / q**2
    d10 = -d11
    return np.array([[d00, d01], [d10, d11]])


def bernoulli_pmf(x, p):
    """Return pmf of Bernoulli distribution evaluated at x."""
    return p**x * (1 - p)**(1 - x)


def bernoulli_pmf_prime(x, p):
    """Return derivative of Bernoulli pmf evaluated at x."""
    return bernoulli_pmf(x, p) * (x/p - (1-x)/(1-p))


# Utility functions
def norm_params(t_dists, e_dists):
    """Return parameters as their normalized values."""
    t_dists_norm = {}
    for s1, t_dist in t_dists.items():
        z_sum = sum([exp(z) for z in t_dist.values()])
        t_dists_norm[s1] = {s2: exp(z)/z_sum for s2, z in t_dist.items()}
    e_dists_norm = {}
    for s, (zp, zpi, zq0, zq1) in e_dists.items():
        e_dists_norm[s] = 1 / (1 + exp(-zp)), 1 / (1 + exp(-zpi)), exp(zq0), exp(zq1)
    return t_dists_norm, e_dists_norm


def get_expectations(t_dists_norm, e_dists_norm, start_dist, record):
    """Return record updated with expected values of states and transitions given model parameters."""
    # Instantiate model
    tree, state_seq, emit_seq = record['tree'], record['state_seq'], record['emit_seq']
    tree_probabilities = {}
    tree_primes_pi, tree_primes_q0, tree_primes_q1 = {}, {}, {}
    for state, (p, pi, q0, q1) in e_dists_norm.items():
        tree_probabilities[state] = utils.get_tree_probability(tree, pi, q0, q1)
        tree_primes_pi[state] = get_tree_prime_pi(tree, q0, q1)
        tree_primes_q0[state] = get_tree_prime_q0(tree, pi, q0, q1)
        tree_primes_q1[state] = get_tree_prime_q1(tree, pi, q0, q1)
    record.update({'tree_probabilities': tree_probabilities,
                   'tree_primes_pi': tree_primes_pi, 'tree_primes_q0': tree_primes_q0, 'tree_primes_q1': tree_primes_q1})
    model = hmm.HMM(t_dists_norm,
                    {state: utils.BinomialArrayRV(p, tree_probabilities[state]) for state, (p, _, _, _) in e_dists_norm.items()},
                     start_dist)

    # Get expectations
    fs, ss_f = model.forward(emit_seq)
    bs, ss_b = model.backward(emit_seq)
    record['nis'] = model.forward_backward1(emit_seq, fs, ss_f, bs, ss_b)
    record['nijs'] = model.forward_backward2(emit_seq, fs, ss_f, bs, ss_b)

    # Calculate likelihood
    px = reduce(lambda x, y: x + y, map(log, ss_f))
    pxy = model.joint_likelihood(emit_seq, state_seq)
    record['ll'] = pxy - px

    return record


eta = 1E-5  # Learning rate
epsilon = 1E-2  # Convergence criterion
iter_num = 200  # Max number of iterations
num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])

if __name__ == '__main__':
    # Load regions
    OGid2regions = {}
    state_set = set()
    with open('labels.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            OGid, start, stop, state = line.rstrip('\n').split('\t')
            if state != '0':  # Skip terminal insertions as actual state
                state_set.add(state)
            try:
                OGid2regions[OGid].append((int(start), int(stop), state))
            except KeyError:
                OGid2regions[OGid] = [(int(start), int(stop), state)]

    # Convert MSAs to records containing state-emissions sequences and other data
    records = []
    for OGid, regions in OGid2regions.items():
        # Load MSA and trim terminal insertions
        msa = read_fasta(f'../../ortho_MSA/realign_hmmer/out/{OGid}.afa')
        msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]
        if regions[-1][2] == '0':
            start, _, _ = regions[-1]
            regions = regions[:-1]
            trim = []
            for header, seq in msa:
                trim.append((header, seq[:start]))
            msa = trim
        if regions[0][2] == '0':
            _, stop, _ = regions[0]
            trim = []
            for header, seq in msa:
                trim.append((header, seq[stop:]))
            msa = trim

            offset = regions[0][1]
            rs = []
            for start, stop, state in regions[1:]:
                rs.append((start-offset, stop-offset, state))
            regions = rs

        # Create emission sequence
        col0 = []
        emit_seq = []
        for j in range(len(msa[0][1])):
            col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
            emit0 = all([c0 == c for c0, c in zip(col0, col)])
            emit_seq.append((emit0, j))  # The tree probabilities are pre-calculated, so emission value is its index
            col0 = col

        # Create state sequence
        state_seq = []
        for (start, stop, state) in regions:
            state_seq.extend((stop - start) * [state])

        # Load tree and convert to vectors at tips
        tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
        tree = tree.shear([spid for spid, _ in msa])
        tips = {tip.name: tip for tip in tree.tips()}
        for spid, seq in msa:
            tip = tips[spid]
            conditional = np.zeros((2, len(seq)))
            for j, sym in enumerate(seq):
                if sym in ['-', '.']:
                    conditional[0, j] = 1
                else:
                    conditional[1, j] = 1
            tip.conditional = conditional

        # Create count dictionaries
        mis = hmm.count_states(state_seq, state_set)
        mijs = hmm.count_transitions(state_seq, state_set)

        records.append({'OGid': OGid, 'tree': tree, 'state_seq': state_seq, 'emit_seq': emit_seq,
                        'mis': mis, 'mijs': mijs})

    # Initialize parameters
    t_dists = {'1A': {'1A': 0.997, '1B': 0.001, '2': 0.001, '3': 0.001},
               '1B': {'1A': 0.001, '1B': 0.997, '2': 0.001, '3': 0.001},
               '2': {'1A': 0.001, '1B': 0.001, '2': 0.997, '3': 0.001},
               '3': {'1A': 0.001, '1B': 0.001, '2': 0.001, '3': 0.997}}
    e_dists = {'1A': (0.95, 0.9, 4, 1),
               '1B': (0.85, 0.75, 2, 1),
               '2': (0.90, 0.5, 3, 3),
               '3': (0.99, 0.01, 0.5, 10)}
    start_dist = {'1A': 0.25, '1B': 0.20, '2': 0.30, '3': 0.25}

    t_dists = {s1: {s2: log(v) for s2, v in t_dist.items()} for s1, t_dist in t_dists.items()}
    e_dists = {s: (-log(1/p-1), -log(1/pi-1), log(q0), log(q1)) for s, (p, pi, q0, q1) in e_dists.items()}

    # Gradient descent
    j, ll0 = 0, None
    models = []
    while j < iter_num:
        # Calculate expectations and likelihoods
        t_dists_norm, e_dists_norm = norm_params(t_dists, e_dists)
        with mp.Pool(processes=num_processes) as pool:
            records = pool.starmap(get_expectations, [(t_dists_norm, e_dists_norm, start_dist, record) for record in records])

        # Save and report updated parameters
        ll = sum([record['ll'] for record in records])
        models.append((ll, t_dists_norm, e_dists_norm))

        print('ITERATION', j)
        print('ll:', ll)
        print('t_dists_norm:', t_dists_norm)
        print('e_dists_norm:', e_dists_norm)
        print()

        # Check convergence
        if ll0 is not None and abs(ll - ll0) < epsilon:
            break
        ll0 = ll

        # Update parameters
        for record in records:
            # Unpack variables
            emit_seq = record['emit_seq']
            mis, mijs = record['mis'], record['mijs']
            nis, nijs = record['nis'], record['nijs']
            tree_probabilities = record['tree_probabilities']
            tree_primes_pi, tree_primes_q0, tree_primes_q1 = record['tree_primes_pi'], record['tree_primes_q0'], record['tree_primes_q1']

            # Update t_dists
            for s1, t_dist in t_dists.items():
                mn_sum = sum([mijs[(s1, s2)] - nijs[(s1, s2)] for s2 in state_set])
                z_sum = sum([exp(z) for z in t_dist.values()])
                for s2, z in t_dist.items():
                    d = -(mijs[(s1, s2)] - nijs[(s1, s2)] - exp(z)/z_sum * mn_sum)  # Equation 2.17
                    t_dist[s2] -= eta * d

            # Update e_dists
            for s, (zp, zpi, zq0, zq1) in e_dists.items():
                p, pi, q0, q1 = 1 / (1 + exp(-zp)), 1 / (1 + exp(-zpi)), exp(zq0), exp(zq1)
                dzp, dzpi, dzq0, dzq1 = 0, 0, 0, 0
                for i, emit in enumerate(emit_seq):
                    # Equations 2.15 and 2.16 (emission parameter phi only)
                    mn = mis[s][i] - nis[s][i]
                    dzp -= mn / bernoulli_pmf(emit[0], p) * bernoulli_pmf_prime(emit[0], p) * p / (1 + exp(zp))
                    dzpi -= mn / tree_probabilities[s][emit[1]] * tree_primes_pi[s][emit[1]] * pi / (1 + exp(zpi))
                    dzq0 -= mn / tree_probabilities[s][emit[1]] * tree_primes_q0[s][emit[1]] * q0
                    dzq1 -= mn / tree_probabilities[s][emit[1]] * tree_primes_q1[s][emit[1]] * q1
                e_dists[s] = (zp - eta * dzp, zpi - eta * dzpi, zq0 - eta * dzq0, zq1 - eta * dzq1)

        j += 1

    # Save parameters
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open('out/model.json', 'w') as file:
        _, t_dists, e_dists = max(models)
        json.dump({'t_dists': t_dists, 'e_dists': e_dists, 'start_dist': start_dist}, file)

"""
NOTES
This HMM uses a two-state phylo-CTMC emission distribution on the gap patterns. It also uses a Bernoulli distribution on
if the pattern of gaps is the same as in the previous column. The parameters are trained discriminatively using gradient
descent.

The gradients are calculated using the formulas in:
Krogh A, Riis SK. Hidden Neural Networks. Neural Computation. 11, 541-563. 1999.

DEPENDENCIES
../../ortho_MSA/realign_hmmer/realign_hmmer.py
    ../../ortho_MSA/realign_hmmer/out/*.afa
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../config/labels.tsv
../insertion_hmm/fit.py
    ../insertion_hmm/out/model.json
"""