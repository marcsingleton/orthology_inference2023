"""Calculate hidden tip values for missing segments."""

import json
import os
import re

import numpy as np
import skbio
from numpy import exp, log
from src.utils import read_fasta
from src.ortho_MSA.phylo import get_transition_matrix


def get_tree_pmf(tree, tip, pi, q0, q1, p0, p1, param):
    """Return probability of tree given tips."""
    s, conditional = get_conditional(tree, tip, q0, q1, p0, p1, param)
    l = ((exp(s) * conditional) * [[1-pi], [pi]]).sum(axis=0)
    return l


def get_conditional(tree, tip, q0, q1, p0, p1, param, inplace=False):
    """Return conditional probabilities of tree given tips and node state."""
    if not inplace:
        tree = tree.copy()  # Make copy so computations do not change original tree

    for node in tree.postorder():
        if node.is_tip():
            node.s = np.zeros(node.value.shape[1])
            if node.name != tip.name:
                m = [[1-p0, p0], [p1, 1-p1]]
            elif param == 'p0':
                m = [[1-p0, p0], [0, 0]]
            elif param == 'p1':
                m = [[0, 0], [p1, 1-p1]]
            else:
                raise ValueError('"param" is not "p0" or "p1"')
            node.conditional = np.matmul(m, node.value)
        else:
            ss, ps = [], []
            for child in node.children:
                s, conditional = child.s, child.conditional
                m = get_transition_matrix(q0, q1, child.length)
                p = np.matmul(m, conditional)

                ss.append(s)
                ps.append(p)

            conditional = np.product(np.stack(ps), axis=0)
            s = conditional.sum(axis=0)
            node.conditional = conditional / s  # Normalize to 1 to prevent underflow
            node.s = log(s) + np.sum(np.stack(ss), axis=0)  # Pass forward scaling constant in log space

    return tree.s, tree.conditional


def get_tip_pmf(tree, tip, pi, q0, q1, p0, p1):
    pmf0 = get_tree_pmf(tree, tip, pi, q0, q1, p0, p1, 'p0')
    pmf1 = get_tree_pmf(tree, tip, pi, q0, q1, p0, p1, 'p1')
    pmf = pmf0 + pmf1
    return np.stack([pmf0, pmf1]) / pmf


ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
gnid_regex = r'gnid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'

tree_template = skbio.read('../../ortho_tree/consensus_GTR2/out/NI.nwk', 'newick', skbio.TreeNode)
cutoff = 0.8

with open('../missing_hmm/out/model.json') as file:
    model_json = json.load(file)
e_dist = model_json['e_dists']['2']
pi, q0, q1, p0, p1 = [e_dist[param] for param in ['pi', 'q0', 'q1', 'p0', 'p1']]

OGids = [path.removesuffix('.tsv') for path in os.listdir('../missing_trim/out/trims/') if path.endswith('.tsv')]
for OGid in OGids:
    # Get missing segments
    ppid2missing = {}
    with open(f'../missing_trim/out/trims/{OGid}.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            missing = []
            for s in fields['slices'].split(','):
                if s:
                    start, stop = s.split('-')
                    missing.append((int(start), int(stop)))
            ppid2missing[fields['ppid']] = missing

    if not any([missing for missing in ppid2missing.values()]):
        continue

    # Load MSA
    msa = []
    for header, seq in read_fasta(f'../insertion_trim/out/trims/{OGid}.afa'):
        ppid = re.search(ppid_regex, header).group(1)
        spid = re.search(spid_regex, header).group(1)
        gnid = re.search(gnid_regex, header).group(1)
        msa.append({'ppid': ppid, 'spid': spid, 'seq': seq})

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

    # Get hidden tip values
    records = []
    for record in msa:
        ppid, spid = record['ppid'], record['spid']
        pmf = get_tip_pmf(tree, tips[spid], pi, q0, q1, p0, p1)
        binary = pmf[1] >= cutoff
        for start, stop in ppid2missing[ppid]:
            records.append({'ppid': ppid, 'gnid': gnid, 'spid': spid,
                            'start': start, 'stop': stop,
                            'seq': ''.join(['1' if b else '0' for b in binary[start:stop]])})

    # Write to file
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with open(f'out/{OGid}.fa', 'w') as file:
        for record in records:
            ppid, gnid, spid, start, stop, seq = [record[key] for key in ['ppid', 'gnid', 'spid', 'start', 'stop', 'seq']]
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(record['seq']), 80)])
            file.write(f'>ppid={ppid}|gnid={gnid}|spid={spid}|start={start}|stop={stop}\n{seqstring}\n')
