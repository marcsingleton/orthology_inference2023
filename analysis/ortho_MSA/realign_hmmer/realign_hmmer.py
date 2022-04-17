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
    run(f'../../../bin/hmmbuild --hand --eset {eset_scalar*gnidnum} --wnone out/{OGid}.hmm ../realign_trim/out/{OGid}.sto > out/{OGid}.txt', shell=True, check=True)
    run(f'../../../bin/hmmalign --outformat afa out/{OGid}.hmm {path} > out/{OGid}_temp.afa', shell=True, check=True)

    # Remove excess gaps
    msa = read_fasta(f'out/{OGid}_temp.afa')
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
    with open(f'out/{OGid}.afa', 'w') as file:
        for header, seq1 in msa:
            seq2 = ''.join([seq1[s] for s in slices])
            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)])
            file.write(f'{header}\n{seqstring}\n')
    os.remove(f'out/{OGid}_temp.afa')


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
eset_scalar = 1.5  # Effective sequence number scalar; multiplies the gnidnum by this value to add weight to observed sequences

OGid2data = {}
with open('../OG_filter/out/OG_filter.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        fields = line.rstrip('\n').split('\t')
        OGid, ppidnum, gnidnum = fields[1], int(fields[3]), int(fields[5])
        OGid2data[OGid] = (ppidnum, gnidnum)

if __name__ == '__main__':
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        pool.map(hmm_align, OGid2data)

"""
NOTES
This version up-weights the number of sequences by 1.5x to prevent spurious alignments. These are typically caused by
unusual patterns of indels and chance similarities between segments that are over-weighted by the strong priors that
HMMer uses when building its profiles. Placing more weight on the observed sequences typically corrects these issues. 

DEPENDENCIES
../make_fastas1/make_fastas1.py
    ../make_fastas1/out/*.fa
../make_fastas2/make_fastas2.py
    ../make_fastas2/out/*.fa
../realign_trim/realign_trim.py
    ../realign_trim/out/*.afa
"""