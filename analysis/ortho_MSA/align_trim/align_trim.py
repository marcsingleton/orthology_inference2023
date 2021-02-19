"""Trim alignments, removing spurious sequences and regions."""

import os
import re
from itertools import chain

import skbio
from trim import trim_msa

# TAXONOMIC CONDITIONS
conditions = [(set(['dinn', 'dgri']), 1),
              (set(['dnov', 'dvir', 'dhyd', 'dmoj', 'dnav']), 4),
              (set(['dper', 'dpse']), 1),
              (set(['dobs', 'dgua', 'dsub']), 2),
              (set(['dana', 'dbip', 'dkik', 'dser']), 3),
              (set(['dele', 'dfic', 'dtak', 'dbia']), 3),
              (set(['deug', 'dere']), 1),
              (set(['dmel']), 1),
              (set(['dmau', 'dsec']), 1)]

if not os.path.exists('out/'):
    os.mkdir('out/')

paths1 = [f'../align_fastas1/out/{file}' for file in os.listdir('../align_fastas1/out/') if file.endswith('.mfa')]
paths2 = [f'../align_fastas2-2/out/{file}' for file in os.listdir('../align_fastas2-2/out/') if file.endswith('.mfa')]
for path in chain(paths1, paths2):
    # Load MSA
    msa1 = skbio.read(path, format='fasta', into=skbio.TabularMSA, constructor=skbio.Protein)
    for seq in msa1:
        match = re.match('ppid=([A-Za-z_0-9]+)\|gnid=([A-Za-z_0-9]+)\|spid=([a-z]+)', seq.metadata['id'])
        d = {key: value for key, value in zip(['ppid', 'gnid', 'spid'], match.groups())}
        seq.metadata.update(d)

    # Test for species groups
    spids = set([seq.metadata['spid'] for seq in msa1])
    if any([len(spids & group) < num for group, num in conditions]):
        continue

    # Trim MSA and save
    msa2 = trim_msa(msa1)
    skbio.write(msa1, format='fasta', into=f'out/{path[-8:]}', max_width=80)

"""
DEPENDENCIES
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
"""