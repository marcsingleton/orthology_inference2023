"""Convert polypeptide alignments to nucleotide alignments."""

import os
import re


def translate(nt_seq, aa_seq):
    delta = len(nt_seq) // 3 - len(aa_seq)
    if abs(delta) > 1:
        print('File:', file_id)
        print('PPID:', ppid)
        print(f'Error: Lengths of translated genomic ({len(nt_seq) // 3}) '
              f'and protein ({len(aa_seq)}) sequences are incompatible.')
        print()
        return
    num_codon = (len(nt_seq) - 3 * max(0, delta)) // 3  # Exclude last codon only if delta == 1

    tr_syms = []  # Translated symbols
    for i in range(0, num_codon):
        codon = nt_seq[3 * i:3 * i + 3]
        tr = 'X' if 'N' in codon else ttable[codon]
        aa = aa_seq[i]
        if tr != aa:
            if i == 0 and aa == 'M':  # Correct non-standard start codon
                tr_syms.append(aa)
                print('File:', file_id)
                print('PPID:', ppid)
                print('Warning: Non-standard start codon detected. '
                      f'Translated residue {ttable[nt_seq[i:i + 3]]} converted to observed residue {aa}.')
                print()
            elif tr == 'X':  # Correct unknown amino acid
                tr_syms.append(aa)
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Warning: Converted unknown translated residue to observed residue {aa}.')
                print()
            elif tr == '*':
                tr_syms.append(aa)
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Warning: Converted stop codon {num_codon-i-1} residues from C-terminus to observed residue {aa}.')
                print()
            else:  # Mismatch
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Error: Mismatched translated ({tr}) and observed ({aa}) residues.')
                print()
                return
        else:
            tr_syms.append(aa)
    if delta == -1:
        tr_syms.append(aa_seq[-1])
        print('File:', file_id)
        print('PPID:', ppid)
        print(f'Warning: Missing final residue {aa_seq[-1]} appended to translated sequence.')
        print()
    return tr_syms


pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+(\.[0-9]+)?)'}

# Parse parameters
params = []
with open('params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, cds_path = line.split()
        params.append((spid, source, cds_path))

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
for spid, source, cds_path in params:
    with open(cds_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid = re.search(pp_regex[source], line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            ppid2cds[ppid] = ''.join(seqlines)

if not os.path.exists('out/'):
    os.mkdir('out/')

for file_id in filter(lambda x: x.endswith('.mfa'), os.listdir('../align_fastas1/out/')):
    # Load alignments
    aa_aligns = []
    with open('../align_fastas1/out/' + file_id) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line
                ppid = re.search('ppid=([NXYPFBpp0-9_\.]+)\|', header)[1]
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            aa_aligns.append((ppid, header, ''.join(seqlines)))

    # Translate CDS
    nt_aligns = []
    for ppid, header, aa_align in aa_aligns:
        aa_seq = aa_align.replace('-', '')
        nt_seq = ppid2cds[ppid]
        tr_syms = translate(nt_seq, aa_seq)
