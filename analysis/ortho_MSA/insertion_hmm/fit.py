"""Fit HMM to segmented alignments using conditional maximum likelihood via gradient descent."""

import json
import multiprocessing as mp
import os
import re
from functools import reduce

import numpy as np
import skbio
import src.ortho_MSA.hmm as hmm
from src.ortho_MSA import phylo
from numpy import exp, log
from src.utils import read_fasta


def norm_params(t_dists, e_dists):
    """Return parameters as their normalized values."""
    t_dists_norm = {}
    for s1, t_dist in t_dists.items():
        z_sum = sum([exp(z) for z in t_dist.values()])
        t_dists_norm[s1] = {s2: exp(z)/z_sum for s2, z in t_dist.items()}
    e_dists_norm = {}
    for s, e_dist in e_dists.items():
        za, zb, zpi, zq0, zq1, zp0, zp1 = [e_dist[param] for param in ['a', 'b', 'pi', 'q0', 'q1', 'p0', 'p1']]
        e_dists_norm[s] = {'a': exp(za), 'b': exp(zb),
                           'pi': 1 / (1 + exp(-zpi)), 'q0': exp(zq0), 'q1': exp(zq1),
                           'p0': 1 / (1 + exp(-zp0)), 'p1': 1 / (1 + exp(-zp1))}
    return t_dists_norm, e_dists_norm


def unnorm_params(t_dists_norm, e_dists_norm):
    """Return parameters as their unnormalized values for gradient descent."""
    t_dists = {}
    for s1, t_dist in t_dists_norm.items():
        t_dists[s1] = {s2: log(v) for s2, v in t_dist.items()}
    e_dists = {}
    for s, e_dist in e_dists_norm.items():
        a, b, pi, q0, q1, p0, p1 = [e_dist[param] for param in ['a', 'b', 'pi', 'q0', 'q1', 'p0', 'p1']]
        e_dists[s] = {'a': log(a), 'b': log(b),
                      'pi': log(pi / (1 - pi)), 'q0': log(q0), 'q1': log(q1),
                      'p0': log(p0 / (1 - p0)), 'p1': log(p1 / (1 - p1))}
    return t_dists, e_dists


def get_gradients(t_dists_norm, e_dists_norm, start_dist, record):
    """Return record updated with expected values of states and transitions given model parameters."""
    # Unpack record fields
    n = record['n']
    tree, state_seq, emit_seq = record['tree'], record['state_seq'], record['emit_seq']
    mis, mijs = record['mis'], record['mijs']

    # Pre-calculate probabilities for each state as array
    betabinom_pmfs, tree_pmfs = {}, {}
    e_dists_rv = {}
    for s, e_dist in e_dists_norm.items():
        a, b, pi, q0, q1, p0, p1 = [e_dist[param] for param in ['a', 'b', 'pi', 'q0', 'q1', 'p0', 'p1']]
        betabinom_pmf = phylo.get_betabinom_pmf(emit_seq, n, a, b)
        tree_pmf = phylo.get_tree_pmf(tree, pi, q0, q1, p0, p1)

        betabinom_pmfs[s] = betabinom_pmf
        tree_pmfs[s] = tree_pmf
        e_dists_rv[s] = phylo.ArrayRV(betabinom_pmf * tree_pmf)

    # Instantiate model and get expectations
    model = hmm.HMM(t_dists_norm, e_dists_rv, start_dist)
    idx_seq = list(range(len(state_seq)))  # Everything is pre-calculated, so emit_seq is the emit index
    fs, ss_f = model.forward(idx_seq)
    bs, ss_b = model.backward(idx_seq)
    nis = model.forward_backward1(idx_seq, fs, ss_f, bs, ss_b)
    nijs = model.forward_backward2(idx_seq, fs, ss_f, bs, ss_b)

    # Calculate likelihood
    px = reduce(lambda x, y: x + y, map(log, ss_f))
    pxy = model.joint_likelihood(idx_seq, state_seq)
    ll = pxy - px

    # Get t_dists gradients
    t_grads = {}
    for s1, t_dist in t_dists_norm.items():
        t_grad = {}
        mn_sum = sum([mijs[(s1, s2)] - nijs[(s1, s2)] for s2 in t_dist])
        for s2, p in t_dist.items():
            t_grad[s2] = -(mijs[(s1, s2)] - nijs[(s1, s2)] - p * mn_sum)  # Equation 2.20
        t_grads[s1] = t_grad

    # Get e_dists gradients
    e_grads = {}
    for s, e_dist in e_dists_norm.items():
        a, b, pi, q0, q1, p0, p1 = [e_dist[param] for param in ['a', 'b', 'pi', 'q0', 'q1', 'p0', 'p1']]
        betabinom_pmf = betabinom_pmfs[s]
        betabinom_prime_a = phylo.get_betabinom_prime(emit_seq, n, a, b, 'a')
        betabinom_prime_b = phylo.get_betabinom_prime(emit_seq, n, a, b, 'b')
        tree_pmf = tree_pmfs[s]
        tree_prime_pi = phylo.get_tree_prime(tree, pi, q0, q1, p0, p1, 'pi')
        tree_prime_q0 = phylo.get_tree_prime(tree, pi, q0, q1, p0, p1, 'q0')
        tree_prime_q1 = phylo.get_tree_prime(tree, pi, q0, q1, p0, p1, 'q1')
        tree_prime_p0 = phylo.get_tree_prime(tree, pi, q0, q1, p0, p1, 'p0')
        tree_prime_p1 = phylo.get_tree_prime(tree, pi, q0, q1, p0, p1, 'p1')

        # Equations 2.15 and 2.16 (emission parameter phi only)
        e_grad = {}
        mn = np.array([mi - ni for mi, ni in zip(mis[s], nis[s])])
        e_grad['a'] = -mn / betabinom_pmf * betabinom_prime_a * a
        e_grad['b'] = -mn / betabinom_pmf * betabinom_prime_b * b
        e_grad['pi'] = -mn / tree_pmf * tree_prime_pi * pi * (1 - pi)
        e_grad['q0'] = -mn / tree_pmf * tree_prime_q0 * q0
        e_grad['q1'] = -mn / tree_pmf * tree_prime_q1 * q1
        e_grad['p0'] = -mn / tree_pmf * tree_prime_p0 * p0 * (1 - p0)
        e_grad['p1'] = -mn / tree_pmf * tree_prime_p1 * p1 * (1 - p1)
        e_grads[s] = e_grad

    return {'ll': ll, 't_grads': t_grads, 'e_grads': e_grads}


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
spid_regex = r'spid=([a-z]+)'

eta = 0.05  # Learning rate
gamma = 0.85  # Momentum
epsilon = 1E-1  # Convergence criterion
iter_max = 300  # Max number of iterations

state_set = {'1A', '1B', '2', '3'}
start_set = {'1A', '1B', '2', '3'}
t_sets = {s1: {s2 for s2 in state_set} for s1 in state_set}
tree_template = skbio.read('../../ortho_tree/consensus_GTR2/out/NI.nwk', 'newick', skbio.TreeNode)

t_pseudo = 0.1  # t_dist pseudocounts
start_pseudo = 0.1  # start_dist pseudocounts
e_dists_initial = {'1A': {'a': 0.9, 'b': 0.1, 'pi': 0.95, 'q0': 0.01, 'q1': 0.01, 'p0': 0.01, 'p1': 0.05},
                   '1B': {'a': 1, 'b': 0.4, 'pi': 0.5, 'q0': 0.2, 'q1': 0.25, 'p0': 0.015, 'p1': 0.05},
                   '2': {'a': 0.7, 'b': 0.15, 'pi': 0.75, 'q0': 0.05, 'q1': 0.075, 'p0': 0.1, 'p1': 0.4},
                   '3': {'a': 1.6, 'b': 0.35, 'pi': 0.01, 'q0': 0.01, 'q1': 0.01, 'p0': 0.025, 'p1': 0.01}}

if __name__ == '__main__':
    # Load labels
    OGid2labels = {}
    label_set = set()
    with open('labels.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            OGid, start, stop, label = fields['OGid'], int(fields['start']), int(fields['stop']), fields['label']
            label_set.add(label)
            try:
                OGid2labels[OGid].append((start, stop, label))
            except KeyError:
                OGid2labels[OGid] = [(start, stop, label)]

    if state_set != label_set:
        raise RuntimeError('label_set is not equal to state_set')

    # Check label validity
    for OGid, labels in OGid2labels.items():
        start0, stop0, label0 = labels[0]
        if start0 != 0:
            print(f'First interval for {OGid} does not begin at 0')
        for start, stop, label in labels[1:]:
            if label0 == label:
                print(f'State for interval ({OGid}, {start}, {stop}) equals previous state')
            if stop0 != start:
                print(f'Start for interval ({OGid}, {start}, {stop}) does not equal previous stop')
            if start >= stop:
                print(f'Start for interval ({OGid}, {start}, {stop}) is greater than stop')
            stop0, label0 = stop, label

    # Convert MSAs to records containing state-emissions sequences and other data
    records = []
    for OGid, labels in OGid2labels.items():
        # Load MSA
        msa = []
        for header, seq in read_fasta(f'../realign_fastas/out/{OGid}.afa'):
            spid = re.search(spid_regex, header).group(1)
            msa.append({'spid': spid, 'seq': seq})

        # Create emission sequence
        column0 = []
        emit_seq = []
        for j in range(len(msa[0]['seq'])):
            column = [1 if msa[i]['seq'][j] in ['-', '.'] else 0 for i in range(len(msa))]
            emit0 = sum([c0 == c for c0, c in zip(column0, column)])
            emit_seq.append(emit0)  # The tree probabilities are pre-calculated, so emission value is its index
            column0 = column
        emit_seq = np.array(emit_seq)

        # Load tree and convert to vectors at tips
        tree = tree_template.shear([record['spid'] for record in msa])
        tips = {tip.name: tip for tip in tree.tips()}
        for record in msa:
            spid, seq = record['spid'], record['seq']
            value = np.zeros((2, len(seq)))
            for j, sym in enumerate(seq):
                if sym in ['-', '.']:
                    value[0, j] = 1
                else:
                    value[1, j] = 1
            tip = tips[spid]
            tip.value = value

        # Create state sequence
        state_seq = []
        for start, stop, label in labels:
            state_seq.extend((stop - start) * [label])

        # Create indicator and count dictionaries
        mis = hmm.state_indicators(state_seq, state_set)
        mijs = hmm.transition_counts(state_seq, t_sets)

        records.append({'OGid': OGid, 'n': len(msa), 'tree': tree, 'state_seq': state_seq, 'emit_seq': emit_seq,
                        'mis': mis, 'mijs': mijs})
    records = sorted(records, key=lambda x: len(x['state_seq']), reverse=True)  # Process the longest alignments first

    # Calculate start_dist from background distribution of states
    state_counts = {s: start_pseudo for s in start_set}
    for labels in OGid2labels.values():
        for start, stop, label in labels:
            state_counts[label] += stop - start
    state_sum = sum(state_counts.values())
    start_dist = {s: count / state_sum for s, count in state_counts.items()}

    # Initialize t_dist from observed transitions
    t_counts = {s1: {s2: t_pseudo for s2 in t_set} for s1, t_set in t_sets.items()}
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

    # Initialize e_dists from initial values
    e_dists_norm = e_dists_initial.copy()

    # Gradient descent
    t_dists, e_dists = unnorm_params(t_dists_norm, e_dists_norm)
    t_momenta = {s1: {s2: None for s2 in t_dist} for s1, t_dist in t_dists.items()}
    e_momenta = {s: {param: None for param in e_dist} for s, e_dist in e_dists.items()}
    ll0 = None
    history = []
    for iter_num in range(1, iter_max + 1):
        # Calculate gradients
        t_dists_norm, e_dists_norm = norm_params(t_dists, e_dists)
        with mp.Pool(processes=num_processes) as pool:
            gradients = pool.starmap(get_gradients, [(t_dists_norm, e_dists_norm, start_dist, record) for record in records])

        # Save and report parameters from previous update
        ll = sum([gradient['ll'] for gradient in gradients])
        history.append({'iter_num': iter_num, 'll': ll, 't_dists_norm': t_dists_norm, 'e_dists_norm': e_dists_norm})

        print(f'ITERATION {iter_num} / {iter_max}')
        print('\tll:', ll)
        print('\tt_dists_norm:', t_dists_norm)
        print('\te_dists_norm:', e_dists_norm)

        # Check convergence
        if iter_num > 1 and abs(ll - ll0) < epsilon:
            break
        ll0 = ll

        # Accumulate and apply gradients
        for s1, t_dist in t_dists.items():
            for s2 in t_dist:
                grad_stack = np.hstack([gradient['t_grads'][s1][s2] for gradient in gradients])
                if iter_num > 1:
                    dz = gamma * t_momenta[s1][s2] + eta * grad_stack.sum() / len(grad_stack)
                else:
                    dz = eta * grad_stack.sum() / len(grad_stack)
                t_dist[s2] -= dz
                t_momenta[s1][s2] = dz

        for s, e_dist in e_dists.items():
            for param in e_dist:
                grad_stack = np.hstack([gradient['e_grads'][s][param] for gradient in gradients])
                if iter_num > 1:
                    dz = gamma * e_momenta[s][param] + eta * grad_stack.sum() / len(grad_stack)
                else:
                    dz = eta * grad_stack.sum() / len(grad_stack)
                e_dists[s][param] -= dz
                e_momenta[s][param] = dz

    # Save history and best model parameters
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open('out/history.json', 'w') as file:
        json.dump(history, file, indent='\t')
    with open('out/model.json', 'w') as file:
        model = max(history, key=lambda x: x['ll'])
        json.dump({'t_dists': model['t_dists_norm'], 'e_dists': model['e_dists_norm'], 'start_dist': start_dist}, file, indent='\t')

"""
NOTES
This HMM uses a two-state phylo-CTMC emission distribution on the gap patterns. It also uses a Bernoulli distribution on
if the pattern of gaps is the same as in the previous column. The parameters are trained discriminatively using gradient
descent.

The gradients are calculated using the formulas in:
Krogh A, Riis SK. Hidden Neural Networks. Neural Computation. 11, 541-563. 1999.

DEPENDENCIES
../../ortho_tree/consensus_GTR2/consensus_GTR2.py
    ../../ortho_tree/consensus_GTR2/out/NI.nwk
../realign_fastas/realign_fastas.py
    ../realign_fastas/out/*.afa
./labels.tsv
"""