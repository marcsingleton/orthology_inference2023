"""Fit HMM to segmented alignments using conditional maximum likelihood via gradient descent."""

import json
import multiprocessing as mp
import os
import re
from functools import reduce

import src.ortho_MSA.hmm as hmm
from numpy import exp, log
from src.utils import read_fasta


class msaBernoulli:
    def __init__(self, ps):
        self.ps = ps

    def pmf(self, x):
        p = self.ps[x[0]]
        if x[1] == 1:
            return p
        else:
            return 1 - p

    def rvs(self, random_state=None):
        # Required for HMM since it has a simulate method
        # Simulations aren't used here, so it's an empty method
        pass


# Gradient functions
def bernoulli_pmf(x, p):
    """Return pmf of Bernoulli distribution evaluated at x."""
    return p**x * (1 - p)**(1 - x)


def bernoulli_pmf_prime(x, p):
    """Return derivative of Bernoulli pmf evaluated at x."""
    return bernoulli_pmf(x, p) * (x/p - (1-x)/(1-p))


# Utility functions
def norm_params(t_dists, e_param):
    """Return parameters to their normalized values."""
    t_dists_norm = {}
    for s1, t_dist in t_dists.items():
        z_sum = sum([exp(z) for z in t_dist.values()])
        t_dists_norm[s1] = {s2: exp(z)/z_sum for s2, z in t_dist.items()}
    e_param_norm = 1 / (1 + exp(-e_param))
    return t_dists_norm, e_param_norm


def get_expectations(t_dists_norm, e_param_norm, start_dist, record):
    """Return record updated with expected values of states and transitions given model parameters."""
    # Instantiate model
    p_seq, state_seq, emit_seq = record['p_seq'], record['state_seq'], record['emit_seq']
    model = hmm.HMM(t_dists_norm,
                    {'0': msaBernoulli(p_seq), '1': msaBernoulli([e_param_norm for _ in range(len(p_seq))])},
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


eta = 1E-2
epsilon = 1E-2
iter_num = 200
ppid_regex = r'ppid=([A-Za-z0-9_]+)'
num_processes = 2

if __name__ == '__main__':
    # Load regions
    OGid2regions = {}
    state_set = set()
    with open('labels.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            OGid, ppid, start, stop, state = fields['OGid'], fields['ppid'], int(fields['start']), int(fields['stop']), fields['state']
            state_set.add(state)
            try:
                OGid2regions[(OGid, ppid)].append((start, stop, state))
            except KeyError:
                OGid2regions[(OGid, ppid)] = [(start, stop, state)]

    # Convert MSAs to records containing state-emissions sequences and other data
    records = []
    for (OGid, ppid), regions in OGid2regions.items():
        # Load MSA and extract seq
        msa = read_fasta(f'../insertion_trim/out/{OGid}.afa')
        seq = [seq for header, seq in msa if re.search(ppid_regex, header).group(1) == ppid][0]

        # Create Bernoulli sequence
        p_seq = []
        for j in range(len(msa[0][1])):
            col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
            p = sum(col) / len(col)
            p_seq.append(p)

        # Create emission sequence
        emit_seq = []
        for j, sym in enumerate(seq):
            if sym in ['-', '.']:
                emit_seq.append((j, 1))
            else:
                emit_seq.append((j, 0))

        # Create state sequence and counts for initializing parameters
        state_seq = []
        e_count = {0: 0, 1: 0}
        for (start, stop, state) in regions:
            state_seq.extend((stop - start) * [state])
            if state == '1':
                count = len([sym for sym in seq[start:stop] if sym in ['-', '.']])
                e_count[1] += count
                e_count[0] += stop - start - count

        # Create count dictionaries
        mis = hmm.count_states(state_seq, state_set)
        mijs = hmm.count_transitions(state_seq, state_set)

        records.append({'OGid': OGid, 'p_seq': p_seq, 'state_seq': state_seq, 'emit_seq': emit_seq,
                        'mis': mis, 'mijs': mijs,
                        'start_state': regions[0][2], 'e_count': e_count})

    # Initialize parameters (compiling counts from individual records plus pseudocounts)
    t_counts = {state: {s: 1 for s in state_set} for state in state_set}
    e_count = {0: 1, 1: 1}
    start_count = {'0': 1, '1': 1}
    for record in records:
        for (state0, state), count in record['mijs'].items():
            t_counts[state0][state] += count
        for s, count in record['e_count'].items():
            e_count[s] += count
        start_count[record['start_state']] += 1

    t_dists = {}
    for state, t_count in t_counts.items():
        total = sum(t_count.values())
        t_dists[state] = {s: count / total for s, count in t_count.items()}
    e_param = e_count[1] / sum(e_count.values())
    total = sum(start_count.values())
    start_dist = {state: count / total for state, count in start_count.items()}

    t_dists = {s1: {s2: log(v) for s2, v in t_dist.items()} for s1, t_dist in t_dists.items()}
    e_param = -log(1/e_param-1)

    # Gradient descent
    j, ll0 = 0, None
    models = []
    while j < iter_num:
        # Calculate expectations and likelihoods
        t_dists_norm, e_param_norm = norm_params(t_dists, e_param)
        with mp.Pool(processes=num_processes) as pool:
            records = pool.starmap(get_expectations, [(t_dists_norm, e_param_norm, start_dist, record) for record in records])

        # Save and report updated parameters
        ll = sum([record['ll'] for record in records])
        models.append((ll, t_dists_norm, e_param_norm))

        print(f'ITERATION {j} / {iter_num}')
        print('\tll:', ll)
        print('\tt_dists_norm:', t_dists_norm)
        print('\te_param_norm:', e_param_norm)

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

            # Update t_dists
            for s1, t_dist in t_dists.items():
                mn_sum = sum([mijs[(s1, s2)] - nijs[(s1, s2)] for s2 in state_set])
                z_sum = sum([exp(z) for z in t_dist.values()])
                for s2, z in t_dist.items():
                    d = -(mijs[(s1, s2)] - nijs[(s1, s2)] - exp(z)/z_sum * mn_sum)
                    t_dist[s2] -= eta * d

            # Update e_dists
            zp = e_param
            p, dzp = 1 / (1 + exp(-zp)), 0
            for i, emit in enumerate(emit_seq):
                mn = mis['1'][i] - nis['1'][i]
                dzp -= mn / bernoulli_pmf(emit[1], p) * bernoulli_pmf_prime(emit[1], p) * p / (1 + exp(zp))
            e_param = zp - eta * dzp

        j += 1

    # Save parameters
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open('out/model.json', 'w') as file:
        _, t_dists, e_param = max(models)
        json.dump({'t_dists': t_dists, 'e_param': e_param, 'start_dist': start_dist}, file)

"""
NOTES
The gradients are calculated using the formulas in:
Krogh A, Riis SK. Hidden Neural Networks. Neural Computation. 11, 541-563. 1999.

DEPENDENCIES
../insertion_trim/extract.py
    ../insertion_trim/out/*.afa
./labels.tsv
"""