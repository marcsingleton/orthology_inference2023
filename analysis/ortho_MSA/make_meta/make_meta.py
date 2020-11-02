"""Make meta alignments from nucleotide alignments of genes."""

import os
from random import randrange

import Bio.AlignIO as AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# Extract column pools
colpools = {'100red1': [], '100red2': [], '50red': [], '0red': []}
for file_id in filter(lambda x: x.endswith('.mfa'), os.listdir('../align_aa2nt1/out/')):
    align = AlignIO.read(f'../align_aa2nt1/out/{file_id}', 'fasta')
    for i in range(len(align[0])):
        col = [(seq.description[-4:], seq[i]) for seq in align]
        gapnum = [seq[i] for seq in align].count('-')
        if gapnum <= 0 * len(align):
            colpools['100red1'].append(col)
        if gapnum <= 0 * len(align) and any([col[0][1] != sym[1] for sym in col]):
            colpools['100red2'].append(col)
        if gapnum <= 0.5 * len(align):
            colpools['50red'].append(col)
        if gapnum <= 1 * len(align):
            colpools['0red'].append(col)

# Make meta alignments
for label, colpool in colpools.items():
    if not os.path.exists(f'out/{label}/'):
        os.makedirs(f'out/{label}/')

    print(f'{label}:', len(colpool))
    for samplenum in range(100):
        sample = [colpool[randrange(len(colpool))] for _ in range(10000)]
        seqs = {}
        for col in sample:
            for spid, sym in col:
                try:
                    seqs[spid].append(sym)
                except KeyError:
                    seqs[spid] = [sym]

        align = MultipleSeqAlignment([SeqRecord(Seq(''.join(seq)), id=spid, description=f'{label}_{samplenum}')
                                      for spid, seq in seqs.items()])
        AlignIO.write(align, f'out/{label}/meta_{samplenum}.fasta', 'fasta')
        AlignIO.write(align, f'out/{label}/meta_{samplenum}.phy', 'phylip-relaxed')

"""
OUTPUT
100red1: 3085154
100red2: 1612236
50red: 3328845
0red: 3780294

DEPENDENCIES
../align_aa2nt1/align_aa2nt1.py
    ../align_aa2nt1/out/*.mfa
"""