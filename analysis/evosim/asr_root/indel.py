"""Calculate ancestral sequence reconstruction at root for indel process."""

import os
from io import StringIO

import numpy as np
import skbio
from scipy.special import gammainc
from scipy.stats import gamma
from src.evosim.asr import get_conditional, get_tree
from src.utils import read_fasta

if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = [path.removesuffix('.iqtree') for path in os.listdir('../asr_indel/out/') if path.endswith('.iqtree')]
for OGid in OGids:
    # Load substitution model
    with open(f'../asr_indel/out/{OGid}.iqtree') as file:
        line = file.readline()
        while line != 'State frequencies: (estimated with maximum likelihood)\n':
            line = file.readline()

        # Load equilibrium frequencies
        freqs = np.zeros(2)
        for _ in range(2):
            line = file.readline()
        for i in range(2):
            freq = float(line.rstrip().split(' = ')[1])
            freqs[i] = freq
            line = file.readline()

        # Load rate matrix
        matrix = np.zeros((2, 2))
        for _ in range(3):
            line = file.readline()
        for i in range(2):
            rates = line.split()
            matrix[i] = [float(rate) for rate in rates[1:]]
            line = file.readline()

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

    # Load rate categories
    # In IQ-TREE, only the shape parameter is fit and the rate parameter beta is set to alpha so the mean of gamma distribution is 1
    # The calculations here directly correspond to equation 10 in Yang. J Mol Evol (1994) 39:306-314.
    # Note the equation has a small typo where the difference in gamma function evaluations should be divided by the probability
    # of that category since technically it is the rate given that category
    with open(f'../asr_indel/out/{OGid}.iqtree') as file:
        line = file.readline()
        while not line.startswith('Model of rate heterogeneity:'):
            line = file.readline()
        num_categories = int(line.rstrip().split(' Gamma with ')[1][0])
        alpha = float(file.readline().rstrip().split(': ')[1])
    igfs = []  # Incomplete gamma function evaluations
    for i in range(num_categories+1):
        x = gamma.ppf(i/num_categories, a=alpha, scale=1/alpha)
        igfs.append(gammainc(alpha+1, alpha*x))
    rates = []
    for i in range(num_categories):
        rate = num_categories * (igfs[i+1] - igfs[i])
        rates.append((rate, 1/num_categories))

    # Load sequence and convert to vectors at tips of tree
    mca = read_fasta(f'../asr_indel/out/{OGid}.afa')
    tips = {tip.name: tip for tip in tree.tips()}
    for header, seq in mca:
        tip = tips[header[1:5]]
        conditional = np.zeros((2, len(seq)))
        for j, sym in enumerate(seq):
            conditional[int(sym), j] = 1
        tip.conditional = conditional

    # Get likelihoods for rate categories
    likelihoods = []
    for rate, prior in rates:
        s, conditional = get_conditional(tree, rate * matrix)
        l = np.expand_dims(freqs, -1) * conditional
        likelihoods.append(np.exp(s) * l * prior)

    likelihoods = np.stack(likelihoods)
    likelihoods = likelihoods / likelihoods.sum(axis=(0, 1))
    np.save(f'out/{OGid}_indel.npy', likelihoods)

"""
NOTES
See notes in aa.py for reasoning for re-calculating rates from alpha.

DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../asr_indel/asr_indel.py
    ../asr_indel/out/*.iqtree
    ../asr_indel/out/*.afa
"""