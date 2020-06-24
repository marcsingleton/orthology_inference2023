"""Parse BLAST results as gene graph."""

import json
import os
import re
from itertools import permutations

header = ['length', 'pident', 'nident', 'gaps',
          'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send',
          'evalue', 'bitscore']

def get_BHs(subjects):
    BHs = []
    ppids = set()
    gnid = ppid2gnid[re.search(pp_regex[params[db_species]], subjects[0][1]).group(1)]
    for subject in subjects:
        # Check if ppid added already
        BH_ppid = re.search(pp_regex[params[db_species]], subject[1]).group(1)
        if BH_ppid in ppids:
            continue
        ppids.add(BH_ppid)

        # Check if gnid has changed
        BH_gnid = ppid2gnid[BH_ppid]
        if BH_gnid != gnid:
            break

        BH = {'BH_ppid': BH_ppid, 'BH_gnid': BH_gnid,
              **{key: val for key, val in zip(header, subject[2:])}}
        BHs.append(BH)

    return BHs


def add_BH(ggraph, query_ppid, query_gnid, BH_ppid, BH_gnid, **kwargs):
    try:
        ggraph[query_gnid][BH_gnid][query_ppid][BH_ppid] = kwargs
    except KeyError as err:
        if err.args[0] == query_gnid:
            ggraph[query_gnid] = {BH_gnid: {query_ppid: {BH_ppid: kwargs}}}
        elif err.args[0] == BH_gnid:
            ggraph[query_gnid][BH_gnid] = {query_ppid: {BH_ppid: kwargs}}
        elif err.args[0] == query_ppid:
            ggraph[query_gnid][BH_gnid][query_ppid] = {BH_ppid: kwargs}


pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'(XP_[0-9]+(\.[0-9]+)?)',
            'YO': r'(YOtr[A-Z]{2}[0-9]+\|orf[0-9]+)'}

# Load pp metadata
ppid2gnid = {}
with open('../ppid2meta/out/ppid2meta.tsv') as infile:
    for line in infile:
        ppid, gnid, _ = line.split()
        ppid2gnid[ppid] = gnid

# Parse parameters
params = {}
with open('params.tsv') as infile:
    fields = infile.readline().split()  # Skip header
    for line in infile:
        species, _, source = line.split()
        params[species] = source

# Parse BLAST results
ggraph = {}
for query_species, db_species in permutations(params.keys(), 2):
    with open(f'../blast_AAA/out/{query_species}/{db_species}.blast') as file:
        query_ppid, subjects = None, []
        line = file.readline()
        while line:
            # Record query
            while line.startswith('#'):
                if line == '# BLASTP 2.10.0+\n' and query_ppid is not None:  # Only add if previous search returned no hits
                    add_BH(ggraph, query_ppid, query_gnid, db_species, None)
                elif line.startswith('# Query:'):
                    query_ppid = re.search(pp_regex[params[query_species]], line).group(1)
                    query_gnid = ppid2gnid[query_ppid]
                line = file.readline()

            # Record hits
            while line and not line.startswith('#'):
                subjects.append(line.split())
                line = file.readline()

            # Add best from hit list
            BHs = get_BHs(sorted(subjects, key=lambda x: float(x[-2]))) if subjects \
                  else [{'BH_ppid': db_species, 'BH_gnid': None}]  # In case last search in file returned no hits
            for BH in BHs:
                add_BH(ggraph, query_ppid, query_gnid, **BH)
            query_ppid, subjects = None, []  # Signals current search was successfully recorded

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write graph to file
with open('out/ggraph.json', 'w') as outfile:
    json.dump(ggraph, outfile, indent=1)

# Write graph as adjacency list to file
with open('out/ggraph.tsv', 'w') as outfile:
    for query_gnid, adj_gns in ggraph.items():
        outfile.write(query_gnid + '\t' + ','.join(adj_gns.keys()) + '\n')

"""
DEPENDENCIES
../blast_AAA/blast_AAA.py
    ../blast_AAA/out/*
../ppid2meta/ppid2meta.py
    ../ppid2meta/out/ppid2meta.tsv
./params.tsv
"""