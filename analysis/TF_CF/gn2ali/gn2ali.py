"""Map Drosophila genes to alignments."""

path_members = '../../../data/EggNOGv5/drosophilidae/7214_members.tsv'
path_lut = '../../../data/FB_lookup_tables/fbgn_fbtr_fbpp_fb_2019_03.tsv'

# Create dictionary relating each polypeptide sequence to an alignment
pp2ali = {}
with open(path_members) as file_members:
    for line in file_members:
        fields = line.rstrip().split('\t')
        alignment = fields[1]
        proteins = fields[4].split(',')
        for protein in proteins:
            if protein[:4] == '7227':
                pp2ali[protein[5:]] = alignment

# Loop through gene, transcript, polypeptide lookup table and use pp2ali dictionary to find corresponding alignments
with open('gn2ali.tsv', 'w') as file_out:
    file_out.write('FlyBase_FBgn\tFlyBase_FBtr\tFlyBase_FBpp\tEggNOGv5_ID\n')  # Write header to output file

    with open(path_lut) as file_lut:
        for line in file_lut:
            if line.startswith('##'):  # Skip headers in input file
                continue

            fields = line.rstrip('\n').split('\t')
            if fields[2] != '':
                try:
                    fields.append(pp2ali[fields[2]])
                except:
                    fields.append('')
            else:
                fields.append('')

            file_out.write('\t'.join(fields) + '\n')

"""
NOTES
Because each polypeptide is represented in only one alignment, the mapping does not need to account for duplicates (see
../EggNOGv5_validation/unique_seq).

DEPENDENCIES
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
../../../data/FB_lookup_tables/
    ../../../data/FB_lookup_tables/fbgn_fbtr_fbpp_fb_2019_03.tsv
"""