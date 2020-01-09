"""Filter alignments on species and count criteria and save filtered lists."""

import os

inpaths = ['../filter_unknown_realign/out/7214_noX_members.tsv', '../../../data/EggNOGv5/drosophilidae/7214_members.tsv']
outpaths = ['Dmel_members.tsv', '10_10_members.tsv', 'equal_+5_members.tsv']

for inpath in inpaths:
    # Create subdirectory to save plots and move to that directory
    root, _ = os.path.splitext(os.path.basename(inpath))  # Get name of member file without extension
    if not os.path.exists('out/' + root):
        os.makedirs('out/' + root)
    os.chdir('out/' + root)

    # Opening and closing files is slow, so minimize by performing once at beginning and end
    outfiles = []
    for out in outpaths:
        outfiles.append(open(out, 'w'))  # 'w' option overwrites existing files

    with open('../../' + inpath) as infile:  # Add additional parent directory relative path since current directory is output subdirectory
        num_filters = [0 for _ in outpaths]
        num_alignments = 0
        for line in infile:
            fields = line.rstrip().split('\t')
            conditions = ['7227' in fields[5].split(','),  # Contains D. melanogaster
                          fields[2] == '10' and fields[3] == '10',  # Contains 10 total sequences and 10 unique species
                          fields[2] == fields[3] and int(fields[3]) > 5]  # Total sequences and unique species are equal and contains greater than 5 species
            num_alignments += 1
            for i, out in enumerate(outfiles):
                if conditions[i]:
                    num_filters[i] += 1
                    out.write(line)

    for out in outfiles:
        out.close()

    print(root.upper())
    print('number of total alignments:', num_alignments)
    for i, out in enumerate(outpaths):
        print()
        print(out)
        print('number of alignments:', num_filters[i])
        print('fraction of alignments: ', num_filters[i] / num_alignments)
    print()

    os.chdir('../..')  # Return to initial directory

"""
OUTPUT
7214_NOX_MEMBERS
number of total alignments: 13686

Dmel_members.tsv
number of alignments: 11239
fraction of alignments:  0.8212041502265088

10_10_members.tsv
number of alignments: 4992
fraction of alignments:  0.3647523016220956

equal_+5_members.tsv
number of alignments: 7865
fraction of alignments:  0.5746748502118953

7214_MEMBERS
number of total alignments: 13686

Dmel_members.tsv
number of alignments: 11245
fraction of alignments:  0.8216425544351893

10_10_members.tsv
number of alignments: 5069
fraction of alignments:  0.3703784889668274

equal_+5_members.tsv
number of alignments: 7800
fraction of alignments:  0.5699254712845243

DEPENDENCIES
../filter_unknown_realign/filter_unknown_realign.py
    ../filter_unknown_realign/out/7214_noX_members.tsv
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
"""