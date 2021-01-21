"""Draw test alignments."""

import os
import re

import pandas as pd
import skbio
from src.draw import draw_alignment


def load_alignment(path):
    MSA = []
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
            MSA.append((spid, seq))
    return MSA


rows = []
with open('../../ortho_cluster3/clique4+_community/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for gnid in gnids:
            rows.append({'CCid': CCid, 'OGid': OGid, 'gnid': gnid})
OGs = pd.DataFrame(rows)

test_genes = pd.read_table('test_genes.tsv')
OG_meta = pd.read_table('../OG_meta/out/OG_meta.tsv').drop(['CCid', 'edgenum'], axis=1)
tree = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree = tree.shear([tip.name for tip in tree.tips() if tip.name != 'sleb'])

if not os.path.exists('out/'):
    os.mkdir('out/')

df = OGs.merge(test_genes, on='gnid', how='right').merge(OG_meta, on='OGid', how='left')
df.to_csv('out/OGids.tsv', sep='\t', index=False)

for record in df.dropna().itertuples():
    if record.sqidnum == record.gnidnum:
        MSA = load_alignment(f'../align_fastas1/out/{record.OGid}.mfa')
    else:
        MSA = load_alignment(f'../align_fastas2-1/out/{record.OGid}.mfa')

    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    MSA = sorted(MSA, key=lambda x: order[x[0]])  # Re-order sequences
    draw_alignment(MSA, f'out/{record.OGid}.png')

"""
DEPENDENCIES
../../../src/draw.py
../../ortho_cluster3/clique4+_community/clique4+_community2.py
    ../../ortho_cluster3/clique4+_community/out/ggraph2/5clique/gclusters.txt
../ortho_tree/consensus_tree/consensus_tree.py
    ../ortho_tree/consensus_tree/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-1/align_fastas2-0.py
    ../align_fastas2-1/out/*.mfa
../OG_meta/OG_meta.py
    ../OG_meta/out/OG_meta.tsv
./test_genes.tsv
"""