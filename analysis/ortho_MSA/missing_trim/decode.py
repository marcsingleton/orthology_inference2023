"""Decode sequences into posterior probabilities."""

import json
import multiprocessing as mp
import os
import re

import homomorph
import numpy as np
import skbio
from src.ortho_MSA import phylo
from src.utils import read_fasta


def decode(OGid, model_json, tree_template):
    # Load MSA
    msa = []
    for header, seq in read_fasta(f'../insertion_trim/out/trims/{OGid}.afa'):
        msa_ppid = re.search(ppid_regex, header).group(1)
        msa_spid = re.search(spid_regex, header).group(1)
        msa.append({'ppid': msa_ppid, 'spid': msa_spid, 'seq': seq})

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

    records = []
    for record in msa:
        spid, ppid = record['spid'], record['ppid']

        # Instantiate model
        e_dists_rv = {}
        for s, e_dist in model_json['e_dists'].items():
            pinv, alpha, pi, q0, q1, p0, p1 = [e_dist[param] for param in ['pinv', 'alpha', 'pi', 'q0', 'q1', 'p0', 'p1']]
            pmf = phylo.get_tip_pmf(tree, tips[spid], pinv, k, alpha, pi, q0, q1, p0, p1)
            e_dists_rv[s] = phylo.ArrayRV(pmf)
        model = homomorph.HMM(model_json['t_dists'], e_dists_rv, model_json['start_dist'])

        # Decode states and write
        idx_seq = list(range(len(msa[0]['seq'])))  # Everything is pre-calculated, so emit_seq is the emit index
        fbs = model.forward_backward(idx_seq)
        records.append((ppid, fbs))

    state_set = sorted(model_json['t_dists'])
    with open(f'out/posteriors/{OGid}.tsv', 'w') as file:
        file.write('ppid\t' + '\t'.join(state_set) + '\n')
        for ppid, fbs in records:
            for fb in zip(*[fbs[state] for state in state_set]):
                file.write(f'{ppid}\t' + '\t'.join([str(v) for v in fb]) + '\n')


num_processes = 10#int(os.environ['SLURM_CPUS_ON_NODE'])
ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'
k = 4

tree_template = skbio.read('../../ortho_tree/consensus_GTR2/out/NI.nwk', 'newick', skbio.TreeNode)

if __name__ == '__main__':
    # Load OGids
    OGids = []
    with open('../realign_fastas/out/errors.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            OGid, error_flag1, error_flag2 = fields['OGid'], fields['error_flag1'], fields['error_flag2']
            if error_flag1 == 'False' and error_flag2 == 'False':
                OGids.append(OGid)

    with open('../missing_hmm/out/model.json') as file:
        model_json = json.load(file)

    if not os.path.exists('out/posteriors/'):
        os.makedirs('out/posteriors/')

    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(decode, [(OGid, model_json, tree_template) for OGid in OGids])

"""
DEPENDENCIES
../../ortho_tree/consensus_GTR2/consensus_GTR2.py
    ../../ortho_tree/consensus_GTR2/out/NI.nwk
../insertion_trim/trim.py
    ../insertion_trim/out/*.afa
../missing_hmm/fit.py
    ../missing_hmm/out/model.json
../realign_fastas/realign_fastas.py
    ../realign_fastas/out/errors.tsv
"""