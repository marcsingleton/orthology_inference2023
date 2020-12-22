"""Make FASTAs of OGs with one unique sequence per gene."""

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

# Load pp metadata
ppid2gnid = {}
with open('../../ortho_search/ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        ppid, gnid, _ = line.split()
        ppid2gnid[ppid] = gnid

# Load seqs
gnid2seqs = {}
for spid0, source, prot_path in params:
    with open(prot_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid0 = re.search(pp_regex[source], line).group(1)
                gnid = ppid2gnid[ppid0]
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq0 = ''.join(seqlines)

            try:
                for _, _, seq1 in gnid2seqs[gnid]:
                    if seq0 == seq1:
                        break
                else:
                    gnid2seqs[gnid].append((ppid0, spid0, seq0))
            except KeyError:
                gnid2seqs[gnid] = [(ppid0, spid0, seq0)]

# Load OGs and OG metadata
OGs = {}
with open('../clique5+_community/out/5clique/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        OGs[OGid] = gnids

OGs_meta = pd.read_table('../OGid2meta/out/OGid2meta.tsv')

# Write sequences
gn27 = OGs_meta['gnidnum'] == 27
sp27 = OGs_meta['spidnum'] == 27
sq27 = OGs_meta['sqidnum'] == 27
OGids = OGs_meta.loc[gn27 & sp27 & sq27, 'OGid']

if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid in OGids:
    with open(f'out/{OGid}.tfa', 'w') as file:
        for gnid in OGs[OGid]:
            ppid, spid, seq = gnid2seqs[gnid][0]
            seqstring = '\n'.join(seq[i:i+80] for i in range(0, len(seq), 80)) + '\n'
            file.write(f'>ppid={ppid}|gnid={gnid}|spid={spid}\n')
            file.write(seqstring)

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_search/ppid2meta/ppid2meta.py
    ../../ortho_search/ppid2meta/out/ppid2meta.tsv
../clique5+_community/clique5+_community.py
    ../clique5+_community/out/5clique/gclusters.txt
../OGid2meta/OGid2meta.py
    ../OGid2meta/out/OGid2meta.tsv
./params.tsv
"""