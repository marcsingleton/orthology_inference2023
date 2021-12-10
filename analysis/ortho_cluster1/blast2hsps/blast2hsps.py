"""Extract HSPs from BLAST results"""

import os
import re
from itertools import permutations


def get_bhsps(ihsps):
    ihsps = sorted(ihsps, key=lambda x: x['bitscore'], reverse=True)

    # Get cutoff
    cutoff = 0
    gnids = set()
    for ihsp in ihsps:
        if ihsp['bitscore'] == ihsps[0]['bitscore']:
            gnids.add(ihsp['sgnid'])  # Add gnid to allowable list if bitscore is equal to highest
        if ihsp['sgnid'] not in gnids:
            cutoff = ihsp['bitscore']  # Find cutoff separately in case an accepted gnid has the same bitscore as the cutoff
            break

    # Extract best hsps
    bhsps = []
    ppids = set()
    for ihsp in ihsps:
        if ihsp['bitscore'] <= cutoff:
            break  # Stop recording hits once bitscore is lower than cutoff
        if ihsp['sppid'] in ppids:
            continue  # Record hit parameters for only best hit to a polypeptide

        bhsps.append(ihsp)
        ppids.add(ihsp['sppid'])  # "Mark" ppid as added to skip in future

    return bhsps


pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'(XP_[0-9]+(\.[0-9]+)?)',
            'YO': r'(YOtr[A-Z]{2}[0-9]+\|orf[0-9]+)'}
columns = {'qppid': str, 'qgnid': str, 'qspid': str,
           'sppid': str, 'sgnid': str, 'sspid': str,
           'length': int, 'nident': int, 'gaps': int,
           'qlen': int, 'qstart': int, 'qend': int,
           'slen': int, 'sstart': int, 'send': int,
           'evalue': float, 'bitscore': float}

# Load seq metadata
ppid2gnid = {}
with open('../ppid2meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, _, _ = line.split()
        ppid2gnid[ppid] = gnid

# Parse genomes
genomes = {}
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, _ = line.split()
        genomes[spid] = source

# Parse BLAST results
ohsps = []
nulls = []
for query_spid, subject_spid in permutations(genomes, 2):
    with open(f'../blast_AAA/out/{query_spid}/{subject_spid}.blast') as file:
        query_ppid, ihsps = None, []
        line = file.readline()
        while line:
            # Record query
            while line.startswith('#'):
                if line == '# BLASTP 2.10.0+\n' and query_ppid is not None:  # Only add if previous search returned no hits
                    nulls.append({'qppid': query_ppid, 'qgnid': query_gnid, 'qspid': query_spid, 'sspid': subject_spid})
                elif line.startswith('# Query:'):
                    query_ppid = re.search(pp_regex[genomes[query_spid]], line).group(1)
                    query_gnid = ppid2gnid[query_ppid]
                line = file.readline()

            # Record hits
            while line and not line.startswith('#'):
                fields = line.rstrip('\n').split('\t')
                subject_ppid = re.search(pp_regex[genomes[subject_spid]], fields[1]).group(1)
                subject_gnid = ppid2gnid[subject_ppid]
                values = [query_ppid, query_gnid, query_spid,
                          subject_ppid, subject_gnid, subject_spid,
                          *fields[2:],
                          False]
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

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write HSPs and nulls to file
with open(f'out/hsps.tsv', 'w') as file:
    file.write('\t'.join(columns) + '\n')
    for ohsp in ohsps:
        file.write('\t'.join([str(ohsp[column]) for column in columns]) + '\n')
with open(f'out/nulls.tsv', 'w') as file:
    file.write('\t'.join(['qppid', 'qgnid', 'qspid', 'sspid']) + '\n')
    for null in nulls:
        file.write('\t'.join([null[column] for column in ['qppid', 'qgnid', 'qspid', 'sspid']]) + '\n')

"""
NOTES
This script makes a few decisions about what exactly a best hit to both a polypeptide and a gene are. Regarding
polypeptides, since it is possible for there to be multiple hits to different segments of the sequence, the "best" is
taken as the one with the highest bitscore and all others are ignored. Since the quality of the resulting orthologous
group is likely better assessed by aligning the sequences anyway, these secondary hits are not of high importance.

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