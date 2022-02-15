"""Plot example alignments segmented via posterior decoding."""

import json
import re

import matplotlib.pyplot as plt
import numpy as np
import skbio
import src.hmm as hmm
import src.draw as draw
from numpy import exp, log
from src.brownian2.trim import trim_terminals
from src.utils import read_fasta


class BinomialArrayRV:
    """Class that allows HMM class to interface with pre-computed probability array.

    The rvs method is only defined so the HMM recognizes it as a proper
    random variable. Since random variates are not needed in this script, the
    body of the method is left blank.
    """
    def __init__(self, p, array):
        self.p = p
        self.array = array

    def pmf(self, x):
        p1 = self.p if x[0] else 1 - self.p
        p2 = self.array[x[1]]
        return p1 * p2

    def rvs(self, size=None):
        pass


def get_tree_probability(tree, pi, q0, q1):
    """Return probability of tree given tips."""
    s, conditional = get_conditional(tree, q0, q1)
    l = ((exp(s) * conditional) * [[1-pi], [pi]]).sum(axis=0)
    return l


def get_conditional(node, q0, q1):
    """Return conditional probabilities of tree given tips and node state."""
    child1, child2 = node.children

    if child1.is_tip():
        s1, conditional1 = 0, child1.conditional
    else:
        s1, conditional1 = get_conditional(child1, q0, q1)
    if child2.is_tip():
        s2, conditional2 = 0, child2.conditional
    else:
        s2, conditional2 = get_conditional(child2, q0, q1)

    p1 = get_transition_matrix(q0, q1, child1.length)
    p2 = get_transition_matrix(q0, q1, child2.length)
    conditional = np.matmul(p1, conditional1) * np.matmul(p2, conditional2)
    s = conditional.sum(axis=0)
    conditional = conditional / s  # Normalize to 1 to prevent underflow
    s = log(s) + s1 + s2  # Pass forward scaling constant in log space

    return s, conditional


def get_transition_matrix(q0, q1, t):
    """Return transition matrix for two-state CTMC."""
    q = q0 + q1
    p00 = (q1 + q0 * exp(-q*t)) / q
    p01 = 1 - p00
    p11 = (q0 + q1 * exp(-q*t)) / q
    p10 = 1 - p11
    return np.array([[p00, p01], [p10, p11]])


# Load model parameters
with open('out/model.json') as file:
    params = json.load(file)


# Load OGids
OGids = set()
with open('../config/segments.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid = line.split()[0]
        OGids.add(OGid)

# Load tree
tree_template = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_template.tips())}

# Plot alignments
for OGid in OGids:
    # Load msa and trim terminal insertions
    msa = trim_terminals(read_fasta(f'../../ortho_MSA/realign_hmmer1/out/{OGid}.mfa'))
    msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]

    # Create emission sequence
    col0 = []
    emits = []
    for j in range(len(msa[0][1])):
        col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
        emit0 = all([c0 == c for c0, c in zip(col0, col)])
        emit1 = sum(col)
        emits.append((emit0, j))  # The tree probabilities are pre-calculated, so emission value is its index
        col0 = col

    # Load tree and convert to vectors at tips
    tree = tree_template.deepcopy().shear([spid for spid, _ in msa])
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

    # Instantiate model
    e_dists_rv = {}
    for state, (p, pi, q0, q1) in params['e_dists'].items():
        array = get_tree_probability(tree, pi, q0, q1)
        e_dists_rv[state] = BinomialArrayRV(p, array)
    model = hmm.HMM(params['t_dists'], e_dists_rv, params['start_dist'])

    # Decode states and plot
    fbs = model.forward_backward(emits)
    msa = [seq.upper() for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only

    draw.plot_msa_lines(msa, [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}_wide.png', bbox_inches='tight')
    plt.close()

    draw.plot_msa_lines(msa, [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(8, 8))
    plt.savefig(f'out/{OGid}_tall.png', bbox_inches='tight')
    plt.close()

"""
DEPENDENCIES
../../ortho_MSA/realign_hmmer1/realign_hmmer1.py
    ../../ortho_MSA/realign_hmmer1/out/*.mfa
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../config/segments.tsv
./fit.py
    ./out/model.json
"""