"""Parse BLAST results as gene graph."""

import json
import os
import re
from itertools import permutations


def get_BHs(subjects):
    BHs = []
    ppids = set()
    gnids = set()
    cutoff = 0

    for subject in subjects:
        BH_ppid = re.search(pp_regex[params[subject_species]], subject[1]).group(1)
        BH_gnid = ppid2gnid[BH_ppid]
        BH_bitscore = float(subject[-1])

        if BH_bitscore == float(subjects[0][-1]):  # Add gnid to allowable list if bitscore is equal to highest
            gnids.add(BH_gnid)
        if BH_gnid not in gnids:  # Record bitscore once a hit is not within allowable list
            cutoff = BH_bitscore
            break

    for subject in subjects:
        BH_ppid = re.search(pp_regex[params[subject_species]], subject[1]).group(1)
        BH_gnid = ppid2gnid[BH_ppid]
        BH_bitscore = float(subject[-1])

        if BH_bitscore <= cutoff:  # Stop recording hits once bitscore is lower than cutoff
            break
        if BH_ppid in ppids:  # Record hit parameters for only best hit to a polypeptide
            continue

        header = ['length', 'pident', 'nident', 'gaps', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send', 'evalue', 'bitscore']
        BH = {'BH_ppid': BH_ppid, 'BH_gnid': BH_gnid, **{key: val for key, val in zip(header, subject[2:])}}
        BHs.append(BH)
        ppids.add(BH_ppid)  # "Mark" ppid as added to skip in future

    if not BHs:  # In case last search in file returned no hits or there are no unique top gene hits
        BHs.append({'BH_ppid': subject_species, 'BH_gnid': 'null'})

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
with open('../ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        ppid, gnid, _ = line.split()
        ppid2gnid[ppid] = gnid

# Parse parameters
params = {}
with open('params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        species, _, source = line.split()
        params[species] = source

# Parse BLAST results
ggraph = {}
for query_species, subject_species in permutations(params.keys(), 2):
    with open(f'../blast_AAA/out/{query_species}/{subject_species}.blast') as file:
        query_ppid, subjects = None, []
        line = file.readline()
        while line:
            # Record query
            while line.startswith('#'):
                if line == '# BLASTP 2.10.0+\n' and query_ppid is not None:  # Only add if previous search returned no hits
                    add_BH(ggraph, query_ppid, query_gnid, subject_species, 'null')
                elif line.startswith('# Query:'):
                    query_ppid = re.search(pp_regex[params[query_species]], line).group(1)
                    query_gnid = ppid2gnid[query_ppid]
                line = file.readline()

            # Record hits
            while line and not line.startswith('#'):
                subjects.append(line.split())
                line = file.readline()

            # Add best from hit list
            BHs = get_BHs(sorted(subjects, key=lambda x: float(x[-1]), reverse=True))
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
    for query_gnid, BH_gnids in ggraph.items():
        outfile.write(query_gnid + '\t' + ','.join(BH_gnids.keys()) + '\n')

"""
NOTES
This script makes a few decisions about what exactly a best hit to both a polypeptide and a gene are. Regarding
polypeptides, since it is possible for there to be multiple hits to different segments of the sequence, the "best" is
taken as the one with the lowest E-value and all others are ignored. Since the quality of the resulting orthologous
group is likely better assessed by aligning the sequences anyway, these secondary hits are not of high importance.
Regarding genes, a strict adherence to the top hit criterion would allow only polypeptide per gene. This can create
mutually exclusive sets of polypeptide hits within a single orthologous group, a more flexible strategy is allowing
multiple best hits within a gene. The most natural extension of the polypeptide criterion is choosing top hits ranked by
evalue until the parent gene of the polypeptide changes. Subsequent hits to the top gene are considered ambiguous and
ignored (similar to how subsequent hits to a polypeptide are ignored in the polypeptide case). This criterion will
exclude some polypeptides that are associated with the best gene, but since these are ranked lower than hits to
polypeptides in other genes, they should not be included since other stronger hits are discarded.

DEPENDENCIES
../blast_AAA/blast_AAA.py
    ../blast_AAA/out/*
../ppid2meta/ppid2meta.py
    ../ppid2meta/out/ppid2meta.tsv
./params.tsv
"""