"""Calculate ancestral sequence reconstruction at root for indel process."""

import os
from io import StringIO

import numpy as np
import scipy.linalg as linalg
import skbio


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