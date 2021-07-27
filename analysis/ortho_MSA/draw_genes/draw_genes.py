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
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, _, sqid = line.split()
        ppid2meta[ppid] = (gnid, sqid)

# Load OGs
rows = []
with open('../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for ppid in ppids:
            gnid, sqid = ppid2meta[ppid]
            rows.append({'CCid': CCid, 'OGid': OGid, 'ppid': ppid, 'gnid': gnid, 'sqid': sqid})
OGs = pd.DataFrame(rows)

# Load pOG metadata and test genes
OG_meta = pd.read_table('../OG_meta/out/OG_meta.tsv')
genes = pd.read_table('genes.tsv')

# Load tree
tree = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree = tree.shear([tip.name for tip in tree.tips() if tip.name != 'sleb'])
order = {tip.name: i for i, tip in enumerate(tree.tips())}

# Draw alignments
if not os.path.exists('out/'):
    os.mkdir('out/')

df = OGs[['gnid', 'OGid']].drop_duplicates().merge(OG_meta, on='OGid', how='right').merge(genes, on='gnid', how='right')
df.to_csv('out/OGs.tsv', sep='\t', index=False)

for row in df.dropna().itertuples():
    if row.sqidnum == row.gnidnum:
        msa = load_msa(f'../align_fastas1/out/{row.OGid}.mfa')
    else:
        msa = load_msa(f'../align_fastas2-2/out/{row.OGid}.mfa')

    msa = [seq[1] for seq in sorted(msa, key=lambda x: order[x[0]])]  # Re-order sequences and extract seq only
    im = draw_msa(msa)
    plt.imsave(f'out/{row.OGid}.png', im)

"""
DEPENDENCIES
../../../src/draw.py
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity2.py
    ../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../OG_meta/OG_meta.py
    ../OG_meta/out/OG_meta.tsv
./genes.tsv
"""