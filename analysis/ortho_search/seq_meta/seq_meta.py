"""Extract and count polypeptide IDs, storing associated metadata as tsv."""

import os
import re

pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}
gn_regex = {'FlyBase': r'parent=(FBgn[0-9]+)',
            'NCBI': r'db_xref=GeneID:([0-9]+)'}

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path, tcds_path = line.split()
        genomes.append((spid, source, prot_path, tcds_path))

# Extract and count polypeptide IDs
ppid_counts = {}  # Counts for each PPID to find duplicates
ppid2meta = {}  # PPID to gene and species
gnid2seqs = {}  # GNID to PPIDs with unique sequences
num_headers = 0
for spid, source, prot_path, tcds_path in genomes:
    # Find parent genes in tcds headers
    with open(tcds_path) as file:
        for line in file:
            if line.startswith('>'):
                num_headers += 1
                gn_match = re.search(gn_regex[source], line)
                pp_match = re.search(pp_regex[source], line)
                try:
                    # First group is entire line, second is first match
                    gnid = gn_match.group(1)
                    ppid = pp_match.group(1)

                    ppid2meta[ppid] = (gnid, spid)
                    ppid_counts[ppid] = ppid_counts.get(ppid, 0) + 1
                except AttributeError:
                    print(line)

    # Find representative sequences in prot files
    with open(prot_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid0 = re.search(pp_regex[source], line).group(1)
                gnid, spid = ppid2meta[ppid0]
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq0 = ''.join(seqlines)

            try:
                for ppid1, seq1 in gnid2seqs[gnid]:
                    if seq0 == seq1:
                        ppid2meta[ppid0] = (gnid, spid, ppid1)
                        break
                else:
                    gnid2seqs[gnid].append((ppid0, seq0))
                    ppid2meta[ppid0] = (gnid, spid, ppid0)
            except KeyError:
                gnid2seqs[gnid] = [(ppid0, seq0)]
                ppid2meta[ppid0] = (gnid, spid, ppid0)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write to file
with open('out/seq_meta.tsv', 'w') as outfile:
    for ppid, meta in ppid2meta.items():
        outfile.write(ppid + '\t' + '\t'.join(meta) + '\n')

print('Total headers:', num_headers)
print('Unique IDs:', len(ppid_counts))

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
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../config/genomes.tsv
"""