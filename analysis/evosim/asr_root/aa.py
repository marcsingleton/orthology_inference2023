"""Calculate ancestral sequence reconstruction at root for amino acid process."""

import os
from io import StringIO

import numpy as np
import scipy.linalg as linalg
import skbio
from scipy.special import gammainc
from scipy.stats import gamma


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


def get_conditional(node, matrix):
    child1, child2 = node.children

    if child1.is_tip():
        s1, conditional1 = 0, child1.conditional
    else:
        s1, conditional1 = get_conditional(child1, matrix)
    if child2.is_tip():
        s2, conditional2 = 0, child2.conditional
    else:
        s2, conditional2 = get_conditional(child2, matrix)

    p1 = linalg.expm(matrix * child1.length)
    p2 = linalg.expm(matrix * child2.length)
    conditional = np.matmul(p1, conditional1) * np.matmul(p2, conditional2)
    s = conditional.sum(axis=0)
    conditional = conditional / s  # Normalize to 1 to prevent underflow
    s = np.log(s) + s1 + s2  # Pass forward scaling constant in log space

    return s, conditional


num_submodels = 4
syms = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
sym2idx = {sym: idx for idx, sym in enumerate(syms)}

# Read substitution model
with open('../asr_aa/50red_D.txt') as file:
    # Parse exchangeability matrix
    matrix = np.zeros((len(syms), len(syms)))
    for i in range(19):
        line = file.readline()
        for j, value in enumerate(line.split()):
            matrix[i + 1, j] = float(value)
            matrix[j, i + 1] = float(value)

    # Parse equilibrium frequencies
    for _ in range(2):
        line = file.readline()
    freqs = np.array([float(value) for value in line.split()])
matrix = freqs * matrix
rate = (freqs * matrix.sum(axis=1)).sum()
matrix = matrix / rate  # Normalize average rate to 1
np.fill_diagonal(matrix, -matrix.sum(axis=1))

if not os.path.exists('out/'):
    os.mkdir('out/')

prefixes = [path.rstrip('.iqtree') for path in os.listdir('../asr_aa/out/') if path.endswith('.iqtree')]
for prefix in prefixes:
    # Load trees
    tree1 = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
    with open(f'../asr_aa/out/{prefix}.iqtree') as file:
        line = file.readline()
        while line != 'Tree in newick format:\n':
            line = file.readline()
        for _ in range(2):
            line = file.readline()
    file = StringIO(line)
    tree2 = skbio.read(file, 'newick', skbio.TreeNode)

    # Find split and calculate ratio
    # Since an arbitrarily rooted tree can split only one of the two true root clades, we check for both when
    # traversing the tree. The root is temporarily set as the parent node of the first of the root clades that is found.
    # As the re-rooted tree has three children at the root, this has the effect of splitting that root clade while
    # preserving the other root clade. For simplicity, the unsplit root clade is designated as nodeA and the other as
    # nodeB. Once these labels are defined, the new root is introduced between these nodes with a small amount of "tree
    # "surgery."
    tree1 = tree1.shear([tip.name for tip in tree2.tips()])
    nodeA, nodeB = tree1.children
    tipsA, tipsB = set([tip.name for tip in nodeA.tips()]), set([tip.name for tip in nodeB.tips()])
    for node in tree2.traverse():
        tips = set([tip.name for tip in node.tips()])
        if tips == tipsA:
            nodeA, nodeB = nodeB, nodeA  # Swap labels so nodeA is the preserved root clade
            ratio = nodeA.length / (nodeA.length + nodeB.length)
            break
        elif tips == tipsB:
            ratio = nodeA.length / (nodeA.length + nodeB.length)
            break
    tree2 = tree2.root_at(node)

    # Insert root node
    nodeB1, nodeB2 = None, None
    for child in tree2.children:
        if child.count(tips=True) == nodeA.count(tips=True):
            nodeA = child
        elif nodeB1 is None:
            nodeB1 = child
        else:
            nodeB2 = child
    lengthA = ratio * nodeA.length
    lengthB = (1-ratio) * nodeA.length

    nodeA = skbio.TreeNode(children=nodeA.children, length=lengthA)
    nodeB = skbio.TreeNode(children=[nodeB1, nodeB2], length=lengthB)
    tree = skbio.TreeNode(children=[nodeA, nodeB])

    # Load rate categories
    with open(f'../asr_aa/out/{prefix}.iqtree') as file:
        line = file.readline()
        while not line.startswith('Proportion of invariable sites:'):
            line = file.readline()
        p = float(line.rstrip().split(': ')[1])
        alpha = float(file.readline().rstrip().split(': ')[1])
    igfs = []  # Incomplete gamma function evaluations
    for i in range(num_submodels+1):
        x = gamma.ppf(i/num_submodels, a=alpha, scale=1/alpha)
        igfs.append(gammainc(alpha+1, alpha*x))
    rates = [(0, 0, p)]  # Normalized rates
    for i in range(num_submodels):
        rate = num_submodels/(1-p) * (igfs[i+1] - igfs[i])
        rates.append((i+1, rate, (1-p)/num_submodels))

    # Load sequence and convert to vectors at base of tree
    msa = load_msa(f'../asr_aa/out/{prefix}.mfa')
    tips = {tip.name: tip for tip in tree.tips()}
    for header, seq in msa:
        tip = tips[header[1:5]]
        conditional = np.zeros((len(syms), len(seq)))
        for j, sym in enumerate(seq):
            if sym in sym2idx:
                i = sym2idx[sym]
                conditional[i, j] = 1
            else:
                conditional[:, j] = 1 / len(syms)
        tip.conditional = conditional

    # Get likelihood for invariable submodel
    # (Background probability for symbol if invariant; otherwise 0)
    likelihood = np.zeros((len(syms), len(msa[0][1])))
    for j in range(len(msa[0][1])):
        invariant = True
        sym0 = msa[0][1][j]
        if sym0 not in sym2idx:
            invariant = False
        else:
            for i in range(1, len(msa)):
                sym = msa[i][1][j]
                if sym != sym0:
                    invariant = False
                    break
        if invariant:
            idx = sym2idx[sym0]
            likelihood[idx, j] = freqs[idx]
    likelihood = likelihood * rates[0][2]  # Multiply by prior for submodel

    # Get likelihoods for rate submodels
    likelihoods = [likelihood]
    for _, rate, prior in rates[1:]:
        s, conditional = get_conditional(tree, rate*matrix)
        l = np.expand_dims(freqs, -1) * conditional
        likelihoods.append(np.exp(s) * l * prior)

    likelihoods = np.stack(likelihoods)
    likelihoods = likelihoods / likelihoods.sum(axis=(0, 1))
    np.save(f'out/{prefix}_aa.npy', likelihoods)

"""
NOTES
The rates for each submodel are calculated manually because some non-invariant submodels are rounded to 0 in IQ-TREE's
output. This results in submodels with zero probability, which introduces problems when normalizing. I felt it was
important to preserve the original model structure when calculating the ASRs, so I decided against merging these
submodels with the invariant submodel.

DEPENDENCIES
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../asr_aa/50red_D.txt
../asr_aa/asr_aa.py
    ../asr_aa/out/*.iqtree
    ../asr_aa/out/*.mfa
"""