"""Make FASTAs of remaining OGs."""

import os
import re

import pandas as pd

pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}

# Parse parameters
params = []
with open('params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        params.append((spid, source, prot_path))

# Load seqs
ppid2seq = {}
for spid0, source, prot_path in params:
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

OGs_meta = pd.read_table('../OGid2meta/out/OGid2meta.tsv')
rclusters = pd.read_table('../reduce_gclusters/out/rclusters.tsv').groupby('OGid')

# Write sequences
gn25 = OGs_meta['gnidnum'] == 25
sp25 = OGs_meta['spidnum'] == 25
sq25 = OGs_meta['sqidnum'] == 25
OGids = OGs_meta.loc[~(gn25 & sp25 & sq25), 'OGid']

if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid in OGids:
    with open(f'out/{OGid}.tfa', 'w') as file:
        for row in rclusters.get_group(OGid).itertuples():
            seq = ppid2seq[row.ppid]
            seqstring = '\n'.join(seq[i:i+80] for i in range(0, len(seq), 80)) + '\n'
            file.write(f'>ppid={row.ppid}|gnid={row.gnid}|spid={row.spid}\n')
            file.write(seqstring)

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../OGid2meta/OGid2meta.py
    ../OGid2meta/out/OGid2meta.tsv
../reduce_gclusters/reduce_gclusters.py
    ../reduce_gclusters/out/rclusters.tsv
./params.tsv
"""