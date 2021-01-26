"""Make FASTAs of OGs with one unique sequence per gene."""

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

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, _ = line.split()
        ppid2meta[ppid] = (gnid, spid)

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

# Load pOGs and pOG metadata
pOGs = {}
with open('../../ortho_cluster3/clique4+_pcommunity/out/5clique/pclusters.txt') as file:
    for line in file:
        _, pOGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        pOGs[pOGid] = ppids
pOG_meta = pd.read_table('../pOG_meta/out/pOG_meta.tsv')

# Write sequences
if not os.path.exists('out/'):
    os.mkdir('out/')

pOGids = pOG_meta.loc[pOG_meta['ppidnum'] == pOG_meta['gnidnum'], 'pOGid']
for pOGid in pOGids:
    with open(f'out/{pOGid}.tfa', 'w') as file:
        for ppid in pOGs[pOGid]:
            seq = ppid2seq[ppid]
            gnid, spid = ppid2meta[ppid]
            seqstring = '\n'.join(seq[i:i+80] for i in range(0, len(seq), 80)) + '\n'
            file.write(f'>ppid={ppid}|gnid={gnid}|spid={spid}\n')
            file.write(seqstring)

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity.py
    ../../ortho_cluster3/clique4+_pcommunity/out/5clique/pclusters.txt
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../config/genomes.tsv
../pOG_meta/pOG_meta.py
    ../pOG_meta/out/pOG_meta.tsv
"""