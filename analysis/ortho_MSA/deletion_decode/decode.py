"""Decode long deletions of trimmed alignments."""

import json
import multiprocessing as mp
import os
import re

import homomorph
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

    def rvs(self, size=None):
        # Required for HMM since it has a simulate method
        # Simulations aren't used here, so it's an empty method
        pass


def decode(OGid, params):
    # Load msa
    msa = read_fasta(f'../insertion_trim/out/{OGid}.afa')

    # Create Bernoulli sequence
    ps = []
    for j in range(len(msa[0][1])):
        col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
        p = sum(col) / len(col)
        ps.append(p)

    # Instantiate model
    e_dists_rv = {'0': msaBernoulli(ps), '1': msaBernoulli([params['e_param'] for _ in range(len(ps))])}
    model = homomorph.HMM(params['t_dists'], e_dists_rv, params['start_dist'])

    # Decode states
    records = []
    for header, seq in msa:
        # Create emission sequence
        emits = []
        for j, sym in enumerate(seq):
            if sym in ['-', '.']:
                emits.append((j, 1))
            else:
                emits.append((j, 0))

        ppid = re.search(ppid_regex, header).group(1)
        fbs = model.forward_backward(emits)
        records.append((ppid, fbs))

    # Write decoded states
    with open(f'out/{OGid}.tsv', 'w') as file:
        file.write('\t'.join(['ppid'] + states) + '\n')
        for ppid, fbs in records:
            for fb in zip(*[fbs[state] for state in states]):
                file.write(ppid + '\t' + '\t'.join([str(v) for v in fb]) + '\n')


ppid_regex = r'ppid=([A-Za-z0-9_]+)'
num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
states = ['0', '1']

if __name__ == '__main__':
    with open('../deletion_hmm/out/model.json') as file:
        params = json.load(file)

    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        args = [(path.removesuffix('.afa'), params) for path in os.listdir('../insertion_trim/out/') if path.endswith('.afa')]
        pool.starmap(decode, args)

"""
DEPENDENCIES
../insertion_trim/extract.py
    ../insertion_trim/out/*.afa
../deletion_hmm/fit.py
    ../deletion_hmm/out/model.json
"""