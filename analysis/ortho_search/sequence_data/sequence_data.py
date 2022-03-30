"""Extract protein and gene IDs from FASTAs, storing data as tsv."""

import os
import re

from src.utils import read_fasta

ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}
gnid_regex = {'FlyBase': r'parent=(FBgn[0-9]+)',
              'NCBI': r'db_xref=GeneID:([0-9]+)'}

# Load genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path, tcds_path = line.rstrip('\n').split('\t')
        genomes.append((spid, source, prot_path, tcds_path))

# Extract and count IDs
ppid_counts = {}
ppid2data = {}
gnid2seqs = {}
for spid, source, prot_path, tcds_path in genomes:
    # Find parent genes in tcds headers
    tcds_fasta = read_fasta(tcds_path)
    for header, _ in tcds_fasta:
        gnid_match = re.search(gnid_regex[source], header)
        ppid_match = re.search(ppid_regex[source], header)
        if gnid_match and ppid_match:
            gnid = gnid_match.group(1)
            ppid = ppid_match.group(1)

            ppid2data[ppid] = (gnid, spid)
            ppid_counts[ppid] = ppid_counts.get(ppid, 0) + 1
        else:
            print(header)

    # Find representative sequences in prot files
    prot_fasta = read_fasta(prot_path)
    for header, seq in prot_fasta:
        ppid_match = re.search(ppid_regex[source], header)
        if ppid_match:
            ppid = ppid_match.group(1)
        else:
            print(header)
            continue

        gnid, spid = ppid2data[ppid]
        if gnid in gnid2seqs:
            seqs = gnid2seqs[gnid]
            if seq in seqs:
                ppid2data[ppid] = (gnid, spid, seqs[seq])
            else:
                seqs[seq] = ppid
                ppid2data[ppid] = (gnid, spid, ppid)
        else:
            gnid2seqs[gnid] = {seq: ppid}
            ppid2data[ppid] = (gnid, spid, ppid)

# Write to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/sequence_data.tsv', 'w') as file:
    file.write('ppid\tgnid\tspid\tsqid\n')
    for ppid, data in ppid2data.items():
        file.write(ppid + '\t' + '\t'.join(data) + '\n')

print('Total headers:', sum(ppid_counts.values()))
print('Unique PPIDs:', len(ppid_counts))

"""
OUTPUT
Total headers: 839414
Unique IDs: 839414

NOTES
dyak has a few unusual entries in its translated_cds file which have no proper NCBI gene ID. These are supplied in the
full database file, so they were simply manually corrected in the source FASTAs.

DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/fasta/dmel-all-translation-r6.38.fasta
../config/genomes.tsv
"""