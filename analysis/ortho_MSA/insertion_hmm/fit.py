"""Fit HMM to segmented alignments using conditional maximum likelihood via gradient descent."""

import json
import multiprocessing as mp
import os
import re
from functools import reduce

import numpy as np
import skbio
import src.ortho_MSA.hmm as hmm
import utils
from numpy import exp, log
from scipy.special import beta, comb, digamma
from src.utils import read_fasta


# Gradient functions
def get_tree_prime_pi(tree, q0, q1, r):
    """Return derivative of probability of tree relative to pi."""
    s, conditional = utils.get_conditional(tree, q0, q1, r)
    d = ((exp(s) * conditional) * [[-1], [1]]).sum(axis=0)  # Broadcasting magic
    return d


def get_tree_prime_q(tree, pi, q0, q1, r, name):
    """Return derivative of probability of tree relative to q."""
    derivative = get_conditional_prime_q(tree, q0, q1, r, name)
    d = (derivative * [[1-pi], [pi]]).sum(axis=0)  # Broadcasting magic
    return d


def get_tree_prime_r(tree, pi, q0, q1, r):
    """Return derivative of probability of tree relative to r."""
    derivative = get_conditional_prime_r(tree, q0, q1, r)
    d = (derivative * [[1-pi], [pi]]).sum(axis=0)  # Broadcasting magic
    return d


def get_conditional_prime_q(node, q0, q1, r, name):
    """Return derivative of conditional probabilities relative to q."""
    # Collect product and derivative of product for each branch
    ps = []
    dps = []
    for child in node.children:
        if child.is_tip():
            conditional = child.conditional
            r_branch = r
            derivative = np.zeros((2, conditional.shape[1]))
        else:
            s, conditional = utils.get_conditional(child, q0, q1, r)
            conditional = exp(s) * conditional  # Un-normalize
            r_branch = 0
            derivative = get_conditional_prime_q(child, q0, q1, r, name)

        m = utils.get_transition_matrix(q0, q1, r_branch, child.length)
        p = np.matmul(m, conditional)

        dm = get_transition_matrix_prime_q(q0, q1, r_branch, child.length, name)
        dp = np.matmul(dm, conditional) + np.matmul(m, derivative)

        ps.append(p)
        dps.append(dp)

    # Assemble products and derivatives into terms and then sum (product rule)
    terms = []
    for i in range(len(node.children)):
        term = [dps[j] if i == j else ps[j] for j in range(len(node.children))]
        term = np.product(np.stack(term), axis=0)
        terms.append(term)
    derivative = np.sum(np.stack(terms), axis=0)

    return derivative


def get_conditional_prime_r(node, q0, q1, r):
    """Return derivative of conditional probabilities relative to r."""
    # Collect product and derivative of product for each branch
    ps = []
    dps = []
    for child in node.children:
        if child.is_tip():
            conditional = child.conditional
            r_branch = r
            derivative = np.zeros((2, conditional.shape[1]))
            dm = get_transition_matrix_prime_r(q0, q1, r_branch, child.length)
        else:
            s, conditional = utils.get_conditional(child, q0, q1, r)
            conditional = exp(s) * conditional  # Un-normalize
            r_branch = 0
            derivative = get_conditional_prime_r(child, q0, q1, r)
            dm = np.array([[0, 0], [0, 0]])

        m = utils.get_transition_matrix(q0, q1, r_branch, child.length)
        p = np.matmul(m, conditional)

        dp = np.matmul(dm, conditional) + np.matmul(m, derivative)

        ps.append(p)
        dps.append(dp)

    # Assemble products and derivatives into terms and then sum (product rule)
    terms = []
    for i in range(len(node.children)):
        term = [dps[j] if i == j else ps[j] for j in range(len(node.children))]
        term = np.product(np.stack(term), axis=0)
        terms.append(term)
    derivative = np.sum(np.stack(terms), axis=0)

    return derivative


def get_transition_matrix_prime_q(q0, q1, r, t, name):
    """Return derivative transition matrix for two-state CTMC relative to q."""
    q = q0 + q1
    t += r
    if name == 'q0':
        d00 = -(q1 + (t * q0 ** 2 + t * q0 * q1 - q1) * exp(-q * t)) / q ** 2
        d01 = -d00
        d11 = q1*(1 - (q0 * t + q1 * t + 1) * exp(-q * t)) / q ** 2
        d10 = -d11
        return np.array([[d00, d01], [d10, d11]])
    elif name == 'q1':
        d00 = q0 * (1 - (t * q0 + t * q1 + 1) * exp(-q * t)) / q ** 2
        d01 = -d00
        d11 = -(q0 + (t * q1 ** 2 + t * q0 * q1 - q0) * exp(-q * t)) / q ** 2
        d10 = -d11
        return np.array([[d00, d01], [d10, d11]])
    else:
        raise ValueError('q is not "q0" or "q1"')


def get_transition_matrix_prime_r(q0, q1, r, t):
    """Return derivative transition matrix for two-state CTMC relative to r."""
    q = q0 + q1
    d00 = -q0 * exp(-q * (t + r))
    d01 = -d00
    d11 = -q1 * exp(-q * (t + r))
    d10 = -d11
    return np.array([[d00, d01], [d10, d11]])


def beta_prime(a, b):
    """Return derivative of beta function relative to its first parameter, a."""
    return beta(a, b) * (digamma(a) - digamma(a + b))


def get_betabinom_prime_a(x, n, a, b):
    """Return derivative of beta-binomial pmf relative to its first shape parameter, a."""
    return comb(n, x) * (beta_prime(x + a, n - x + b) * beta(a, b) - beta(x + a, n - x + b) * beta_prime(a, b)) / (beta(a, b))**2


def get_betabinom_prime_b(x, n, a, b):
    """Return derivative of beta-binomial pmf relative to its second shape parameter, b."""
    return comb(n, x) * (beta_prime(n - x + b, x + a) * beta(a, b) - beta(x + a, n - x + b) * beta_prime(b, a)) / (beta(a, b))**2


# Utility functions
def norm_params(t_dists, e_dists):
    """Return parameters as their normalized values."""
    t_dists_norm = {}
    for s1, t_dist in t_dists.items():
        z_sum = sum([exp(z) for z in t_dist.values()])
        t_dists_norm[s1] = {s2: exp(z)/z_sum for s2, z in t_dist.items()}
    e_dists_norm = {}
    for s, params in e_dists.items():
        za, zb, zpi, zq0, zq1, zr = params
        e_dists_norm[s] = exp(za), exp(zb), 1 / (1 + exp(-zpi)), exp(zq0), exp(zq1), exp(zr)
    return t_dists_norm, e_dists_norm


def unnorm_params(t_dists_norm, e_dists_norm):
    """Return parameters as their unnormalized values for gradient descent."""
    t_dists = {}
    for s1, t_dist in t_dists_norm.items():
        t_dists[s1] = {s2: log(v) for s2, v in t_dist.items()}
    e_dists = {}
    for s, params in e_dists_norm.items():
        a, b, pi, q0, q1, r = params
        e_dists[s] = log(a), log(b), log(pi / (1 - pi)), log(q0), log(q1), log(r)
    return t_dists, e_dists


def get_gradients(t_dists_norm, e_dists_norm, start_dist, state_set, record):
    """Return record updated with expected values of states and transitions given model parameters."""
    # Unpack record fields
    n = record['n']
    tree, state_seq, emit_seq = record['tree'], record['state_seq'], record['emit_seq']
    mis, mijs = record['mis'], record['mijs']

    # Pre-calculate probabilities and derivatives for each state as array
    betabinom_pmfs = {}
    betabinom_primes_a, betabinom_primes_b = {}, {}
    tree_pmfs = {}
    tree_primes_pi, tree_primes_q0, tree_primes_q1, tree_primes_r = {}, {}, {}, {}
    e_dists_rv = {}
    for s, params in e_dists_norm.items():
        a, b, pi, q0, q1, r = params
        betabinom_pmfs[s] = utils.get_betabinom_pmf(emit_seq, n, a, b)
        betabinom_primes_a[s] = get_betabinom_prime_a(emit_seq, n, a, b)
        betabinom_primes_b[s] = get_betabinom_prime_b(emit_seq, n, a, b)
        tree_pmfs[s] = utils.get_tree_pmf(tree, pi, q0, q1, r)
        tree_primes_pi[s] = get_tree_prime_pi(tree, q0, q1, r)
        tree_primes_q0[s] = get_tree_prime_q(tree, pi, q0, q1, r, 'q0')
        tree_primes_q1[s] = get_tree_prime_q(tree, pi, q0, q1, r, 'q1')
        tree_primes_r[s] = get_tree_prime_r(tree, pi, q0, q1, r)
        e_dists_rv[s] = utils.ArrayRV(betabinom_pmfs[s] * tree_pmfs[s])

    # Instantiate model and get expectations
    model = hmm.HMM(t_dists_norm, e_dists_rv, start_dist)
    emit_seq = list(range(len(emit_seq)))  # Everything is pre-calculated, so emit_seq is the emit index
    fs, ss_f = model.forward(emit_seq)
    bs, ss_b = model.backward(emit_seq)
    nis = model.forward_backward1(emit_seq, fs, ss_f, bs, ss_b)
    nijs = model.forward_backward2(emit_seq, fs, ss_f, bs, ss_b)

    # Calculate likelihood
    px = reduce(lambda x, y: x + y, map(log, ss_f))
    pxy = model.joint_likelihood(emit_seq, state_seq)
    ll = pxy - px

    # Get t_dists gradients
    t_grads = {}
    for s1, t_dist in t_dists_norm.items():
        t_grad = {}
        mn_sum = sum([mijs[(s1, s2)] - nijs[(s1, s2)] for s2 in state_set])
        for s2, p in t_dist.items():
            t_grad[s2] = -(mijs[(s1, s2)] - nijs[(s1, s2)] - p * mn_sum)  # Equation 2.20
        t_grads[s1] = t_grad

    # Get e_dists gradients
    e_grads = {}
    for s, params in e_dists_norm.items():
        # Equations 2.15 and 2.16 (emission parameter phi only)
        e_grad = {}
        a, b, pi, q0, q1, r = params
        mn = np.array([mi - ni for mi, ni in zip(mis[s], nis[s])])
        e_grad['za'] = (-mn / betabinom_pmfs[s] * betabinom_primes_a[s] * a).sum()
        e_grad['zb'] = (-mn / betabinom_pmfs[s] * betabinom_primes_b[s] * b).sum()
        e_grad['zpi'] = (-mn / tree_pmfs[s] * tree_primes_pi[s] * pi * (1 - pi)).sum()
        e_grad['zq0'] = (-mn / tree_pmfs[s] * tree_primes_q0[s] * q0).sum()
        e_grad['zq1'] = (-mn / tree_pmfs[s] * tree_primes_q1[s] * q1).sum()
        e_grad['zr'] = (-mn / tree_pmfs[s] * tree_primes_r[s] * r).sum()
        e_grads[s] = e_grad

    return {'ll': ll, 't_grads': t_grads, 'e_grads': e_grads}


eta = 1E-5  # Learning rate
epsilon = 1E-2  # Convergence criterion
iter_num = 50  # Max number of iterations
num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'

if __name__ == '__main__':
    # Load labels
    OGid2labels = {}
    state_set = set()
    with open('labels.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            OGid, start, stop, label = fields['OGid'], int(fields['start']), int(fields['stop']), fields['label']
            state_set.add(label)
            try:
                OGid2labels[OGid].append((start, stop, label))
            except KeyError:
                OGid2labels[OGid] = [(start, stop, label)]
    state_set = sorted(state_set)

    # Convert MSAs to records containing state-emissions sequences and other data
    records = []
    for OGid, labels in OGid2labels.items():
        # Load MSA
        msa = []
        for header, seq in read_fasta(f'../realign_hmmer/out/mafft/{OGid}.afa'):
            spid = re.search(spid_regex, header).group(1)
            msa.append({'spid': spid, 'seq': seq})

        # Create emission sequence
        col0 = []
        emit_seq = []
        for j in range(len(msa[0]['seq'])):
            col = [1 if msa[i]['seq'][j] in ['-', '.'] else 0 for i in range(len(msa))]
            emit0 = sum([c0 == c for c0, c in zip(col0, col)])
            emit_seq.append(emit0)  # The tree probabilities are pre-calculated, so emission value is its index
            col0 = col
        emit_seq = np.array(emit_seq)

        # Load tree and convert to vectors at tips
        tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
        tree = tree.shear([record['spid'] for record in msa])
        tips = {tip.name: tip for tip in tree.tips()}
        for record in msa:
            spid, seq = record['spid'], record['seq']
            conditional = np.zeros((2, len(seq)))
            for j, sym in enumerate(seq):
                if sym in ['-', '.']:
                    conditional[0, j] = 1
                else:
                    conditional[1, j] = 1
            tip = tips[spid]
            tip.conditional = conditional

        # Create state sequences
        state_seq = []
        for start, stop, label in labels:
            state_seq.extend((stop - start) * [label])

        # Create count dictionaries
        mis = hmm.count_states(state_seq, state_set)
        mijs = hmm.count_transitions(state_seq, state_set)

        records.append({'OGid': OGid, 'n': len(msa), 'tree': tree, 'state_seq': state_seq, 'emit_seq': emit_seq,
                        'mis': mis, 'mijs': mijs})

    # Calculate start_dist from background distribution of states
    state_counts = {s: 0.1 for s in state_set}
    for labels in OGid2labels.values():
        for start, stop, label in labels:
            state_counts[label] += stop - start
    state_sum = sum(state_counts.values())
    start_dist = {s: count / state_sum for s, count in state_counts.items()}

    # Initialize t_dist from observed transitions
    t_counts = {s1: {s2: 0.1 for s2 in state_set} for s1 in state_set}
    for labels in OGid2labels.values():
        start, stop, label0 = labels[0]
        t_counts[label0][label0] += stop - start - 1
        for start, stop, label1 in labels:
            t_counts[label0][label1] += 1
            t_counts[label1][label1] += stop - start - 1
            label0 = label1
    t_dists_norm = {}
    for s1, t_count in t_counts.items():
        t_sum = sum(t_count.values())
        t_dists_norm[s1] = {s2: count / t_sum for s2, count in t_count.items()}

    # Initialize e_dists
    e_dists_norm = {'1A': (1, 1, 0.95, 0.1, 0.1, 0.001),
                    '1B': (1, 1, 0.75, 1, 1, 0.001),
                    '2': (1, 1, 0.5, 0.5, 0.5, 0.001),
                    '3': (1, 1, 0.01, 0.1, 0.1, 0.001)}

    t_dists, e_dists = unnorm_params(t_dists_norm, e_dists_norm)

    # Gradient descent
    j, ll0 = 0, None
    history = []
    while j <= iter_num:
        # Calculate expectations and likelihoods
        t_dists_norm, e_dists_norm = norm_params(t_dists, e_dists)
        with mp.Pool(processes=num_processes) as pool:
            gradients = pool.starmap(get_gradients, [(t_dists_norm, e_dists_norm, start_dist, state_set, record) for record in records])

        # Save and report parameters from previous update
        ll = sum([gradient['ll'] for gradient in gradients])
        history.append({'iter_num': j, 'll': ll, 't_dists_norm': t_dists_norm, 'e_dists_norm': e_dists_norm})

        print(f'ITERATION {j} / {iter_num}')
        print('\tll:', ll)
        print('\tt_dists_norm:', t_dists_norm)
        print('\te_dists_norm:', e_dists_norm)

        # Check convergence
        if ll0 is not None and abs(ll - ll0) < epsilon:
            break
        ll0 = ll

        # Accumulate and apply gradients
        for s1, t_dist in t_dists.items():
            for s2 in t_dist:
                d = sum([gradient['t_grads'][s1][s2] for gradient in gradients])
                t_dist[s2] -= eta * d

        for s, params in e_dists.items():
            za, zb, zpi, zq0, zq1, zr = params
            dza = sum([gradient['e_grads'][s]['za'] for gradient in gradients])
            dzb = sum([gradient['e_grads'][s]['zb'] for gradient in gradients])
            dzpi = sum([gradient['e_grads'][s]['zpi'] for gradient in gradients])
            dzq0 = sum([gradient['e_grads'][s]['zq0'] for gradient in gradients])
            dzq1 = sum([gradient['e_grads'][s]['zq1'] for gradient in gradients])
            dzr = sum([gradient['e_grads'][s]['zr'] for gradient in gradients])
            e_dists[s] = (za - eta * dza,
                          zb - eta * dzb,
                          zpi - eta * dzpi,
                          zq0 - eta * dzq0,
                          zq1 - eta * dzq1,
                          zr - eta * dzr)

        j += 1

    # Save history and best model parameters
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open('out/history.json', 'w') as file:
        json.dump(history, file)

    with open('out/model.json', 'w') as file:
        model = max(history, key=lambda x: x['ll'])
        json.dump({'t_dists': model['t_dists_norm'], 'e_dists': model['e_dists_norm'], 'start_dist': start_dist}, file)

"""
NOTES
This HMM uses a two-state phylo-CTMC emission distribution on the gap patterns. It also uses a Bernoulli distribution on
if the pattern of gaps is the same as in the previous column. The parameters are trained discriminatively using gradient
descent.

The gradients are calculated using the formulas in:
Krogh A, Riis SK. Hidden Neural Networks. Neural Computation. 11, 541-563. 1999.

DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../insertion_hmm/fit.py
    ../insertion_hmm/out/model.json
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/mafft/*.afa
./labels.tsv
"""