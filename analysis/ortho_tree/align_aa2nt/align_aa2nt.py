"""Convert polypeptide alignments to nucleotide alignments."""

import os
import re
import sys


def get_codons(nt_seq, aa_seq):
    delta = len(nt_seq) // 3 - len(aa_seq)
    if abs(delta) > 1:
        print('File:', file_id)
        print('PPID:', ppid)
        print(f'Error: Lengths of translated genomic ({len(nt_seq) // 3}) '
              f'and protein ({len(aa_seq)}) sequences are incompatible.')
        print()
        return
    num_codon = (len(nt_seq) - 3 * max(0, delta)) // 3  # Exclude last codon only if delta == 1

    codons = []  # Translated symbols
    for i in range(0, num_codon):
        codon = nt_seq[3*i:3*i+3]
        tr = 'X' if 'N' in codon else ttable[codon]
        aa = aa_seq[i]
        if tr != aa:
            if i == 0 and aa == 'M':  # Correct non-standard start codon
                codons.append(codon)
                print('File:', file_id)
                print('PPID:', ppid)
                print('Warning: Non-standard start codon detected. '
                      f'Translated residue M{i+1} converted to observed residue {aa}.')
                print()
            elif tr == 'X':  # Correct unknown amino acid
                codons.append(codon)
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Warning: Converted unknown translated residue X{i+1} to observed residue {aa}.')
                print()
            elif tr == '*':
                codons.append(codon)
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Warning: Converted stop codon *{i+1} to observed residue {aa}.')
                print()
            else:  # Mismatch
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Error: Mismatched translated residue {tr}{i+1} and observed residue {aa}.')
                print()
                return
        else:
            codons.append(codon)
    if delta == -1:
        codons.append(nt_seq[3*(i+1):] + (3-len(nt_seq) % 3) * '-')
        print('File:', file_id)
        print('PPID:', ppid)
        print(f'Warning: Missing final residue {aa_seq[-1]} appended to translated sequence.')
        print()
    return codons


ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, _, cds_path = line.split()
        genomes.append((spid, source, cds_path))

# Load translation table
ttable = {}
with open('ttable.txt') as file:
    lines = [line.rstrip().split(' = ')[1] for line in file]
    for i in range(len(lines[0])):
        aa = lines[0][i]
        codon = ''.join([lines[j][i] for j in range(2, 5)])
        ttable[codon] = aa

# Load CDSs
ppid2cds = {}
for spid, source, cds_path in genomes:
    with open(cds_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid = re.search(ppid_regex[source], line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            ppid2cds[ppid] = ''.join(seqlines)

if not os.path.exists('out/'):
    os.mkdir('out/')

sys.stdout = open('out/out.txt', 'w')  # Redirect stdout to file
for file_id in filter(lambda x: x.endswith('.mfa'), os.listdir('../align_fastas/out/')):
    # Load alignments
    aa_aligns = []
    with open('../align_fastas/out/' + file_id) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            aa_aligns.append((header, ''.join(seqlines)))

    # Translate and write CDS
    nt_aligns = []
    for header, aa_align in aa_aligns:
        ppid = re.search(r'ppid=([NXYPFBp0-9_.]+)\|', header)[1]
        aa_seq = aa_align.replace('-', '')
        nt_seq = ppid2cds[ppid]

        codons = get_codons(nt_seq, aa_seq)
        if not codons:
            print(f'Conversion of {file_id} failed due to {ppid}.')
            print()
            break

        nt_align = []
        for aa in aa_align:
            if aa == '-':
                nt_align.append('---')
            else:
                nt_align.append(codons.pop(0))
        nt_aligns.append((header, ''.join(nt_align)))
    else:
        with open('out/' + file_id, 'w') as file:
            for header, nt_align in nt_aligns:
                alignstring = '\n'.join(nt_align[i:i+80] for i in range(0, len(nt_align), 80)) + '\n'
                file.write(header + alignstring)
sys.stdout.close()

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_cds_from_genomic.fna
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-CDS-r6.34.fasta
../config/genomes.tsv
../align_fastas/align_fastas.py
    ../align_fastas/out/*.mfa
./ttable.txt
"""