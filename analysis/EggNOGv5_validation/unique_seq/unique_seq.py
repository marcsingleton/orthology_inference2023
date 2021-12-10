"""Calculate the distribution of number of alignments for each polypeptide sequence."""

path = '../../../data/EggNOGv5/drosophilidae/7214_members.tsv'

# Create dictionary relating each polypeptide sequence to a list of its alignments
pp2ali = {}
num_seqs = 0
with open(path) as file:
    for line in file:
        fields = line.rstrip('\n').split('\t')
        alignment = fields[1]
        proteins = fields[4].split(',')
        for protein in proteins:
            num_seqs += 1
            try:
                pp2ali[protein].append(alignment)
            except KeyError:
                pp2ali[protein] = [alignment]

# Count number of alignments for each polypeptide
counts = {}
for alis in pp2ali.values():
    count = len(alis)
    counts[count] = counts.get(count, 0) + 1

print('total polypeptide sequences')
print(num_seqs)
print()
print('distribution of number of alignments for each polypeptide sequence')
print(counts)

"""
OUTPUT
total polypeptide sequences
121579

distribution of number of alignments for each polypeptide sequence
{1: 121579}

NOTES
Every polypeptide sequence is represented in only one alignment.

DEPENDENCIES
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
"""