"""Re-align sequences using HMMER."""

import multiprocessing as mp
import os
from subprocess import run


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


def hmm_align(file):
    OGid = file[:-4]
    if os.path.exists(f'../align_fastas1/out/{OGid}.mfa'):
        path = f'../make_fastas1/out/{OGid}.tfa'
    else:
        path = f'../make_fastas2-2/out/{OGid}.tfa'
    run(f'../../../bin/hmmbuild --symfrac 0 --enone --wnone out/{OGid}.hmm ../align_trim/out/{file} > out/{OGid}.txt', shell=True, check=True)
    run(f'../../../bin/hmmalign --outformat afa out/{OGid}.hmm {path}> out/{OGid}_temp.mfa', shell=True, check=True)

    # Remove excess gaps
    msa = load_msa(f'out/{OGid}_temp.mfa')
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
    with open(f'out/{OGid}.mfa', 'w') as file:
        for header, seq1 in msa:
            seq2 = []
            for s in slices:
                seq2.extend(seq1[s])
            seq2 = ''.join(seq2)

            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)]) + '\n'
            file.write(header)
            file.write(seqstring)
    os.remove(f'out/{OGid}_temp.mfa')


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])

if not os.path.exists('out/'):
    os.mkdir('out/')

if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        pool.map(hmm_align, os.listdir('../align_trim/out/'))

"""
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../align_trim/align_trim.py
    ../align_trim/out/*.mfa
"""