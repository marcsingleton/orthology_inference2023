"""Extract and count polypeptide IDs, storing associated data as tsv."""

import os
import re

from src.utils import read_fasta

ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}
gnid_regex = {'FlyBase': r'parent=(FBgn[0-9]+)',
              'NCBI': r'db_xref=GeneID:([0-9]+)'}

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path, tcds_path = line.split()
        genomes.append((spid, source, prot_path, tcds_path))

# Extract and count polypeptide IDs
counts = {}  # Counts for each PPID to find duplicates
ppid2data = {}  # PPID to gene and species
gnid2seqs = {}  # GNID to PPIDs with unique sequences
for spid, source, prot_path, tcds_path in genomes:
    # Find parent genes in tcds headers
    tcds_fasta = read_fasta(tcds_path)
    for header, _ in tcds_fasta:
        gn_match = re.search(gnid_regex[source], header)
        pp_match = re.search(ppid_regex[source], header)
        try:
            # First group is entire line, second is first match
            gnid = gn_match.group(1)
            ppid = pp_match.group(1)
            ppid2data[ppid] = (gnid, spid)
        except AttributeError:
            print(header)

    # Find representative sequences in prot files
    prot_fasta = read_fasta(prot_path)
    for header, seq in prot_fasta:
        ppid = re.search(ppid_regex[source], header).group(1)
        gnid, spid = ppid2data[ppid]
        counts[ppid] = counts.get(ppid, 0) + 1
        if gnid in gnid2seqs:
            seqs = gnid2seqs[gnid]  # Mapping from sequence to PPID
            if seq in seqs:
                ppid2data[ppid] = (gnid, spid, seqs[seq])
            else:
                seqs[seq] = ppid
                ppid2data[ppid] = (gnid, spid, ppid)
        else:
            gnid2seqs[gnid] = {seq: ppid}
            ppid2data[ppid] = (gnid, spid, ppid)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write to file
with open('out/sequence_data.tsv', 'w') as file:
    file.write('ppid\tgnid\tspid\tsqid\n')
    for ppid, data in ppid2data.items():
        file.write(ppid + '\t' + '\t'.join(data) + '\n')

print('Total headers:', sum(counts.values()))
print('Unique IDs:', len(counts))

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