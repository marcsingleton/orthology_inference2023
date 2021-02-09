"""Draw test alignments."""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import skbio
from src.draw import draw_msa


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                spid = re.search(r'spid=([a-z]+)', line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((spid, seq))
    return msa


# Load seq metadata
ppid2gnid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, _, _ = line.split()
        ppid2gnid[ppid] = gnid

# Load pOGs
rows = []
with open('../../ortho_cluster3/subcluster_pgraph/out/pclusters.txt') as file:
    for line in file:
        OGid, pOGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for ppid in ppids:
            rows.append({'OGid': OGid, 'pOGid': pOGid, 'ppid': ppid, 'gnid': ppid2gnid[ppid]})
pOGs = pd.DataFrame(rows)

# Load pOG metadata and test genes
pOG_meta = pd.read_table('../pOG_meta/out/pOG_meta.tsv', dtype={'pCCid': str})
genes = pd.read_table('genes.tsv')

# Load tree
tree = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree = tree.shear([tip.name for tip in tree.tips() if tip.name != 'sleb'])
order = {tip.name: i for i, tip in enumerate(tree.tips())}

# Draw alignments
if not os.path.exists('out/'):
    os.mkdir('out/')

df = pOGs[['gnid', 'OGid', 'pOGid']].drop_duplicates().merge(pOG_meta.drop('OGid', axis=1), on='pOGid', how='right').merge(genes, on='gnid', how='right')
df.to_csv('out/pOGids.tsv', sep='\t', index=False)

for record in df.dropna().itertuples():
    if record.ppidnum == record.gnidnum:
        msa = load_msa(f'../align_fastas1/out/{record.pOGid}.mfa')
    else:
        msa = load_msa(f'../align_fastas2-2/out/{record.pOGid}.mfa')

    msa = [seq[1] for seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa)
    plt.imsave(f'out/{record.pOGid}.png', im)

"""
DEPENDENCIES
../../../src/draw.py
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_cluster3/subcluster_pgraph.py
    ../../ortho_cluster3/subcluster_pgraph/out/pclusters.txt
../ortho_tree/consensus_tree/consensus_tree.py
    ../ortho_tree/consensus_tree/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../pOG_meta/pOG_meta.py
    ../pOG_meta/out/pOG_meta.tsv
./genes.tsv
"""