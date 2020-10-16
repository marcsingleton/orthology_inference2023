"""Align FASTAs of OGs with one unique sequence per gene."""

import multiprocessing as mp
import os
from subprocess import run
from time import time_ns


def get_cmd(file_id):
    cmd = (f'../../../bin/probalign ../make_fastas1/out/{file_id}.tfa '
           f'1> out/{file_id}.mfa 2> out/{file_id}.err')
    return cmd


def run_cmd(file_id, cmd):
    t0 = time_ns()
    run(cmd, shell=True, check=True)
    t1 = time_ns()

    return file_id, str(t1-t0)


num_processes = int(os.environ['SLURM_NTASKS'])

if __name__ == '__main__':
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        cmds = [(file_id, get_cmd(file_id[:-4])) for file_id in os.listdir('../make_fastas1/out/')]
        rows = pool.starmap(run_cmd, cmds)

        # Save results
    with open('out/times.tsv', 'w') as file:
        file.write('file_id\ttime_ns\n')
        for row in rows:
            file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../make_fastas1/make_fastas1.py
    ../make_fastas1/out/*.tfa
"""