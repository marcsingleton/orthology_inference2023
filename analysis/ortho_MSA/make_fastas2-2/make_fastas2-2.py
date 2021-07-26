"""Make FASTAs of remaining OGs using trees."""

import os
import re

import pandas as pd

pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        genomes.append((spid, source, prot_path))

# Load seqs
ppid2seq = {}
for spid, source, prot_path in genomes:
    with open(prot_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid = re.search(pp_regex[source], line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            ppid2seq[ppid] = seq

# Load clusters and pOG metadata
rclusters = pd.read_table('../reduce_tree/out/rclusters.tsv').groupby('OGid')
OG_meta = pd.read_table('../OG_meta/out/OG_meta.tsv')

# Write sequences
if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = OG_meta.loc[~(OG_meta['sqidnum'] == OG_meta['gnidnum']), 'OGid']
for OGid in OGids:
    with open(f'out/{OGid}.tfa', 'w') as file:
        rows = sorted(rclusters.get_group(OGid).itertuples(), key=lambda x: x.spid)
        for row in rows:
            seq = ppid2seq[row.ppid]
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]) + '\n'
            file.write(f'>ppid={row.ppid}|gnid={row.gnid}|spid={row.spid}\n')
            file.write(seqstring)

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../config/genomes.tsv
../OG_meta/OG_meta.py
    ../OG_meta/out/OG_meta.tsv
../reduce_tree/reduce_tree.py
    ../reduce_tree/out/rclusters.tsv
"""