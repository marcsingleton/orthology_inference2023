"""Extract HSPs from BLAST results"""

import multiprocessing as mp
import os
import re
from itertools import groupby, permutations


def parse_blast(query_spid, subject_spid):
    ohsps = []
    nulls = []
    with open(f'../blast_AAA/out/{query_spid}/{subject_spid}.blast') as file:
        query_ppid, ihsps = None, []
        line = file.readline()
        while line:
            # Record query
            while line.startswith('#'):
                if line == '# BLASTP 2.10.1+\n' and query_ppid is not None:  # Only add if previous search returned no hits
                    nulls.append({'qppid': query_ppid, 'qgnid': query_gnid, 'qspid': query_spid, 'sspid': subject_spid})
                elif line.startswith('# Query:'):
                    query_ppid = re.search(pp_regex[genomes[query_spid]], line).group(1)
                    query_gnid = ppid2gnid[query_ppid]
                line = file.readline()

            # Record hits
            while line and not line.startswith('#'):
                fields = line.split()
                subject_ppid = re.search(pp_regex[genomes[subject_spid]], fields[1]).group(1)
                subject_gnid = ppid2gnid[subject_ppid]
                values = [query_ppid, query_gnid, query_spid,
                          subject_ppid, subject_gnid, subject_spid,
                          *fields[2:],
                          False, False]
                ihsps.append({column: f(value) for (column, f), value in zip(columns.items(), values)})
                line = file.readline()

            # Add best from HSP list (and catch cases where the last search returned no HSPs)
            if ihsps:
                ohsps.extend(get_bhsps(ihsps))
            else:
                nulls.append({'qppid': query_ppid, 'qgnid': query_gnid, 'qspid': query_spid, 'sspid': subject_spid})

            if line.startswith('# BLAST processed'):
                break
            query_ppid, ihsps = None, []  # Signals current search was successfully recorded

    # Make output directories
    if not os.path.exists(f'out/hsps/{query_spid}/'):
        os.makedirs(f'out/hsps/{query_spid}/')  # Recursive folder creation
    if not os.path.exists(f'out/nulls/{query_spid}/'):
        os.makedirs(f'out/nulls/{query_spid}/')  # Recursive folder creation

    # Write HSPs and nulls to file
    with open(f'out/hsps/{query_spid}/{subject_spid}.tsv', 'w') as file:
        file.write('\t'.join(columns) + '\n')
        for ohsp in ohsps:
            file.write('\t'.join([str(ohsp[column]) for column in columns]) + '\n')
    with open(f'out/nulls/{query_spid}/{subject_spid}.tsv', 'w') as file:
        file.write('\t'.join(['qppid', 'qgnid', 'qspid', 'sspid']) + '\n')
        for null in nulls:
            file.write('\t'.join([null[column] for column in ['qppid', 'qgnid', 'qspid', 'sspid']]) + '\n')


def get_bhsps(ihsps):
    ihsps = sorted(ihsps, key=lambda x: (x['sppid'], x['sgnid']))
    groups = {key: list(group) for key, group in groupby(ihsps, lambda x: (x['sppid'], x['sgnid']))}
    bs_maxs = sorted([(ppid, gnid, max([ihsp['bitscore'] for ihsp in ihsps])) for (ppid, gnid), ihsps in groups.items()],
                     key=lambda x: x[2], reverse=True)

    # Get cutoff
    cutoff = 0
    gnids = set()
    for ppid, gnid, bs_max in bs_maxs:
        if bs_max == bs_maxs[0][2]:
            gnids.add(gnid)  # Add gnid to allowable list if bitscore is equal to highest
        if gnid not in gnids:
            cutoff = bs_max  # Choose cutoff as the maximum bitscore associated with the next best gene
            break

    # Extract groups
    keys = []
    for ppid, gnid, bs_max in bs_maxs:
        if bs_max <= cutoff:  # Stop recording hits once bitscore is lower than cutoff
            break
        keys.append((ppid, gnid))

    # Optimize HSPs
    bhsps = []
    for key in keys:
        group = sorted(groups[key], key=lambda x: x['bitscore'], reverse=True)
        group[0]['index_hsp'] = True  # Mark "best" HSP as index for subsequent filtering
        chsps = []
        for ihsp in group:
            if not has_overlap(ihsp, chsps, True) and not has_overlap(ihsp, chsps, False):
                ihsp['disjoint'] = True  # Mark HSPs as disjoint greedily, beginning with highest score
                chsps.append(ihsp)
            bhsps.append(ihsp)  # Track all HSPs

    return bhsps


def has_overlap(ihsp, ohsps, query=True):
    if query:
        start_key = 'qstart'
        end_key = 'qend'
    else:
        start_key = 'sstart'
        end_key = 'send'

    for ohsp in ohsps:
        if not (ihsp[start_key] > ohsp[end_key] or ohsp[start_key] > ihsp[end_key]):
            return True
    return False


pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}
columns = {'qppid': str, 'qgnid': str, 'qspid': str,
           'sppid': str, 'sgnid': str, 'sspid': str,
           'length': int, 'nident': int, 'gaps': int,
           'qlen': int, 'qstart': int, 'qend': int,
           'slen': int, 'sstart': int, 'send': int,
           'evalue': float, 'bitscore': float,
           'index_hsp': bool, 'disjoint': bool}
num_processes = int(os.environ['SLURM_NTASKS'])

# Load seq metadata
ppid2gnid = {}
with open('../seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, _, _ = line.split()
        ppid2gnid[ppid] = gnid

# Parse genomes
genomes = {}
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, _, _ = line.split()
        genomes[spid] = source

# Parse BLAST results
if __name__ == '__main__':
    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(parse_blast, permutations(genomes, 2))

"""
NOTES
This script makes a few decisions about what exactly a best hit to both a polypeptide and a gene are. Regarding
polypeptides, since it is possible for there to be multiple hits to different segments of the sequence, the "best" is
taken as the HSP with the largest bitscore. However, all other HSPs from this hit are recorded. Proceeding from highest
to lowest bitscore, these HSPs are marked as disjoint if they do not overlap with any HSPs previously marked as
disjoint. Additionally other polypeptide hits are recorded if the strongest HSP in that hit exceeds the bitscore of
the strongest HSP to a polypeptide that is not associated with the same gene as the best polypeptide hit. In other
words, this approach takes a broad "gene-centric" approach to initially filtering the raw BLAST output where all hits
associated with the best gene are recorded as long as the bitscore of their strongest HSP exceeds that of the next best
gene. If there are multiple HSPs with the same best bitscore, then they are both treated as the best hit. This means a
single polypeptide can have hits to multiple genes which in turn contain multiple polypeptides which in turn contain
multiple HSPs.

Regarding genes, a strict adherence to the top hit criterion would allow only polypeptide per gene. This can create
mutually exclusive sets of polypeptide hits within a single orthologous group, a more flexible strategy is allowing
multiple best hits within a gene. The most natural extension of the polypeptide criterion is choosing top hits ranked by
evalue until the parent gene of the polypeptide changes. Subsequent hits to the top gene are considered ambiguous and
ignored (similar to how subsequent hits to a polypeptide are ignored in the polypeptide case). This criterion will
exclude some polypeptides that are associated with the best gene, but since these are ranked lower than hits to
polypeptides in other genes, they should not be included since other stronger hits are discarded.

DEPENDENCIES
../config/genomes.tsv
../blast_AAA/blast_AAA.py
    ../blast_AAA/out/*
../seq_meta/seq_meta.py
    ../seq_meta/out/seq_meta.tsv
"""