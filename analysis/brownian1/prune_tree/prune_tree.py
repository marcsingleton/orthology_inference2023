""""Prune unneeded lineages from Drosophila tree and save image of pruned tree."""

import Bio.Phylo as Phylo
import matplotlib.pyplot as plt
import os

# Input variables
path = '../../../data/drosophila_tree_25/drosophila-25spec-tree.tre'
remain = set(['ananassae', 'erecta', 'grimshawi',
              'melanogaster', 'mojavensis', 'persimilis',
              'pseudoobscura', 'virilis', 'yakuba', 'willistoni'])

# Read data and filter
tree = Phylo.read(path, 'newick')
remove = set([terminal.name for terminal in tree.get_terminals()]) - remain
for leaf in remove:
    tree.prune(leaf)

# Make output directory and save
if not os.path.exists('out/'):
    os.mkdir('out/')
Phylo.write(tree, 'out/drosophila-10spec-tree.nwk', 'newick')

# Draw tree
fig, ax = plt.subplots(figsize=(7, 4))
Phylo.draw(tree, show_confidence=False, axes=ax, do_show=False)
plt.savefig('out/drosophila-10spec-tree.png')

"""
DEPENDENCIES
../../../data/drosophila_tree_25/drosophila-25spec-tree.tre
"""