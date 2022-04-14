"""Remove duplicate sequences within genes from FASTAs."""

import os
import re

from src.utils import read_fasta

ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+(\.[0-9]+)?)'}

# Load genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path, _ = line.rstrip('\n').split('\t')
        genomes.append((spid, source, prot_path))

# Load sequence data
ppid2sqid = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, _, _, sqid = line.rstrip('\n').split('\t')
        ppid2sqid[ppid] = sqid

if not os.path.exists('out/'):
    os.mkdir('out/')

for spid, source, prot_path in genomes:
    input_fasta = read_fasta(prot_path)
    output_fasta = []

    for header, seq in input_fasta:
        ppid = re.search(ppid_regex[source], header).group(1)
        sqid = ppid2sqid[ppid]
        if ppid == sqid:
            output_fasta.append((header, seq))

    with open(f'out/{spid}.fa', 'w') as file:
        for header, seq in output_fasta:
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
            file.write(f'{header}\n{seqstring}\n')

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.45_FB2022_02/fasta/dmel-all-translation-r6.45.fasta
../config/genomes.tsv
"""