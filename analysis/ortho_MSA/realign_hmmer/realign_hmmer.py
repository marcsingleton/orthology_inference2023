"""Re-align sequences using HMMER."""

import multiprocessing as mp
import os
from subprocess import run

from src.utils import read_fasta


def hmm_align(OGid):
    ppidnum, gnidnum = OGid2data[OGid]
    if ppidnum == gnidnum:
        path = f'../make_fastas1/out/{OGid}.fa'
    else:
        path = f'../make_fastas2/out/{OGid}.fa'
    run(f'../../../bin/hmmbuild --hand --eset {eset_scalar*gnidnum} --wnone out/hmmer/{OGid}.hmm ../cnn_trim/out/{OGid}.sto > out/hmmer/{OGid}.txt', shell=True, check=True)
    run(f'../../../bin/hmmalign --outformat afa out/hmmer/{OGid}.hmm {path} > out/hmmer/{OGid}_temp.afa', shell=True, check=True)

    # Remove excess gaps
    msa = read_fasta(f'out/hmmer/{OGid}_temp.afa')
    slices, idx = [], None
    for j in range(len(msa[0][1])):
        for i in range(len(msa)):
            sym = msa[i][1][j]
            if sym not in ['-', '.']:
                if idx is None:  # Store position only if new slice is not started
                    idx = j
                break
        else:
            if idx is not None:
                slices.append(slice(idx, j))
                idx = None
    if idx is not None:  # Add final slice to end
        slices.append(slice(idx, len(msa[0][1])))

    # Write to file and remove temp alignment
    with open(f'out/hmmer/{OGid}.afa', 'w') as file:
        for header, seq1 in msa:
            seq2 = ''.join([seq1[s] for s in slices])
            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)])
            file.write(f'{header}\n{seqstring}\n')
    os.remove(f'out/hmmer/{OGid}_temp.afa')

    # Load new alignment and identify unaligned regions
    msa1 = trim_terminals(read_fasta(f'out/hmmer/{OGid}.afa'))
    slices, idx = [], None
    for j in range(len(msa1[0][1])):
        for i in range(len(msa1)):
            sym = msa1[i][1][j]
            if sym == '.' or sym.islower():
                if idx is None:  # Store position only if new slice is not started
                    idx = j
                break
        else:
            if idx is not None:
                slices.append(slice(idx, j))
                idx = None
    if idx is not None:  # Add final slice to end
        slices.append(slice(idx, len(msa1[0][1])))

    # Align subsequences and stitch together results
    submsas = []
    for s in slices:
        # Collect subsequences
        subseqs = []
        for header, seq in msa1:
            subseq = seq[s]
            for sym in subseq:
                if sym != '.':
                    subseqs.append((header, subseq.replace('.', '')))
                    break

        # Align subsequences
        if len(subseqs) == 1:
            output = subseqs
        else:
            with open(f'out/mafft/{OGid}_temp.fa', 'w') as file:
                for header, seq in subseqs:
                    seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
                    file.write(f'{header}\n{seqstring}\n')

            cmd = (f'../../../bin/mafft --globalpair --maxiterate 1000 '
                   f'--thread 1 --anysymbol --allowshift --unalignlevel 0.4 --leavegappyregion '
                   f'out/mafft/{OGid}_temp.fa '
                   f'1> out/mafft/{OGid}_temp.afa 2> out/mafft/{OGid}_temp.err')
            run(cmd, shell=True, check=True)

            output = read_fasta(f'out/mafft/{OGid}_temp.afa')
            os.remove(f'out/mafft/{OGid}_temp.fa')
            os.remove(f'out/mafft/{OGid}_temp.afa')
            os.remove(f'out/mafft/{OGid}_temp.err')
        length = len(output[0][1])
        output = {header: seq for header, seq in output}

        # Fill missing subsequences with gaps
        submsa = []
        for header, _ in msa1:
            if header in output:
                submsa.append((header, output[header]))
            else:
                submsa.append((header, length * '-'))
        submsas.append(submsa)

    # Stitch together results
    msa2, idx = [(header, []) for header, _ in msa1], 0
    for s, submsa in zip(slices, submsas):
        for (_, subseq), (_, seq1), (_, seq2) in zip(submsa, msa1, msa2):
            seq2.append(seq1[idx:s.start])
            seq2.append(subseq)
        idx = s.stop
    for (_, seq1), (_, seq2) in zip(msa1, msa2):  # Add final region
        seq2.append(seq1[idx:])

    # Write to file
    with open(f'out/mafft/{OGid}.afa', 'w') as file:
        for header, seq in msa2:
            seq = ''.join(seq)
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
            file.write(f'{header}\n{seqstring}\n')


def trim_terminals(msa):
    """Return MSA with unaligned terminal regions removed.

    Unaligned regions are created by HMMer a symbol in the sequence does not
    match any of the consensus columns columns in the profile HMM. These
    regions are indicated with lowercase symbols for . for gaps.

    Parameters
    ----------
    msa: list of (header, seq)

    Returns
    -------
    msa: list of (header, seq)
    """
    idx = 0
    for j in range(len(msa[0][1])):
        for i in range(len(msa)):
            sym = msa[i][1][j]
            if sym == '.' or sym.islower():
                break
        else:
            idx = j
            break  # if no break exit
    msa = [(header, seq[idx:]) for header, seq in msa]

    idx = len(msa[0][1])
    for j in range(len(msa[0][1]), 0, -1):
        for i in range(len(msa)):
            sym = msa[i][1][j - 1]
            if sym == '.' or sym.islower():
                break
        else:
            idx = j
            break  # if no break exit
    msa = [(header, seq[:idx]) for header, seq in msa]

    return msa


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
eset_scalar = 1.5  # Effective sequence number scalar; multiplies the gnidnum by this value to add weight to observed sequences

OGid2data = {}
with open('../OG_filter/out/OG_filter.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGid, ppidnum, gnidnum = fields['OGid'], int(fields['ppidnum']), int(fields['gnidnum'])
        OGid2data[OGid] = (ppidnum, gnidnum)

if __name__ == '__main__':
    if not os.path.exists('out/hmmer/'):
        os.makedirs('out/hmmer/')
    if not os.path.exists('out/mafft/'):
        os.makedirs('out/mafft/')

    with mp.Pool(processes=num_processes) as pool:
        pool.map(hmm_align, OGid2data)

"""
NOTES
This version up-weights the number of sequences to prevent spurious alignments. These are typically caused by unusual
patterns of indels and chance similarities between segments that are over-weighted by the strong priors that HMMer uses
when building its profiles. Placing more weight on the observed sequences typically corrects these issues. 

DEPENDENCIES
../cnn_trim/cnn_trim.py
    ../cnn_trim/out/*.sto
../make_fastas1/make_fastas1.py
    ../make_fastas1/out/*.fa
../make_fastas2/make_fastas2.py
    ../make_fastas2/out/*.fa
"""