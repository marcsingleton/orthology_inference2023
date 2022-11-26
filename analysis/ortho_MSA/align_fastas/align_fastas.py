"""Align FASTAs of OGs with one unique sequence per gene."""

import multiprocessing as mp
import os
from subprocess import run, CalledProcessError
from time import time_ns


def run_cmd(OGid):
    cmd = (f'../../../bin/mafft --globalpair --maxiterate 1000 '
           f'--thread 1 --anysymbol --allowshift --unalignlevel 0.4 --leavegappyregion '
           f'../make_fastas/out/{OGid}.fa '
           f'1> out/{OGid}.afa 2> out/{OGid}.err')
    try:
        t0 = time_ns()
        run(cmd, shell=True, check=True)
        t1 = time_ns()
        return OGid, str(t1 - t0)
    except CalledProcessError:
        return OGid, 'NaN'


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])

if __name__ == '__main__':
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        OGids = [path.removesuffix('.fa') for path in os.listdir('../make_fastas/out/') if path.endswith('.fa')]
        rows = pool.map(run_cmd, OGids)

    with open('out/times.tsv', 'w') as file:
        file.write('OGid\ttime_ns\n')
        for row in rows:
            file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../make_fastas/make_fastas.py
    ../make_fastas/out/*.fa
"""