"""Calculate ancestral sequence reconstruction at root for indel process."""

import os
from io import StringIO

import numpy as np
import skbio
from asr import get_conditional, get_tree


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


if not os.path.exists('out/'):
    os.mkdir('out/')

prefixes = [path.rstrip('.iqtree') for path in os.listdir('../asr_indel/out/') if path.endswith('.iqtree')]
for prefix in prefixes:
    # Load substitution model
    with open(f'../asr_indel/out/{prefix}.iqtree') as file:
        line = file.readline()
        while line != 'State frequencies: (estimated with maximum likelihood)\n':
            line = file.readline()

        # Parse equilibrium frequencies
        freqs = np.zeros(2)
        for _ in range(2):
            line = file.readline()
        for i in range(2):
            freq = float(line.rstrip().split(' = ')[1])
            freqs[i] = freq
            line = file.readline()

        # Parse rate matrix
        matrix = np.zeros((2, 2))
        for _ in range(3):
            line = file.readline()
        for i in range(2):
            rates = line.split()
            matrix[i] = [float(rate) for rate in rates[1:]]
            line = file.readline()

    # Load trees
    tree1 = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
    with open(f'../asr_indel/out/{prefix}.iqtree') as file:
        line = file.readline()
        while line != 'Tree in newick format:\n':
            line = file.readline()
        for _ in range(2):
            line = file.readline()
    tree2 = skbio.read(StringIO(line), 'newick', skbio.TreeNode)
    tree = get_tree(tree1, tree2)

    # Load rate categories
    rates = []
    with open(f'../asr_indel/out/{prefix}.iqtree') as file:
        line = file.readline()
        while line != ' Category  Relative_rate  Proportion\n':
            line = file.readline()
        for _ in range(4):
            line = file.readline()
            fields = line.split()
            rates.append((fields[0], float(fields[1]), float(fields[2])))

    # Load sequence and convert to vectors at base of tree
    msa = load_msa(f'../asr_indel/out/{prefix}.mfa')
    tips = {tip.name: tip for tip in tree.tips()}
    for header, seq in msa:
        tip = tips[header[1:5]]
        conditional = np.zeros((2, len(seq)))
        for j, sym in enumerate(seq):
            conditional[int(sym), j] = 1
        tip.conditional = conditional

    # Get likelihoods for rate submodels
    likelihoods = []
    for _, rate, prior in rates:
        s, conditional = get_conditional(tree, rate*matrix)
        l = np.expand_dims(freqs, -1) * conditional
        likelihoods.append(np.exp(s) * l * prior)

    likelihoods = np.stack(likelihoods)
    likelihoods = likelihoods / likelihoods.sum(axis=(0, 1))
    np.save(f'out/{prefix}_indel.npy', likelihoods)

"""
DEPENDENCIES
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../asr_indel/asr_indel.py
    ../asr_indel/out/*.iqtree
    ../asr_indel/out/*.mfa
"""