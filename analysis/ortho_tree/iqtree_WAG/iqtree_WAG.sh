# Run IQ-TREE on meta alignments using WAG model

if [ ! -d out/${1}/ ]; then
  mkdir -p out/${1}/
fi

../../../bin/iqtree -s ../make_metaAA/out/${1}/${2}.fasta -m WAG+FO${3}+G -pre out/${1}/${2} -quiet

# DEPENDENCIES
# ../make_metaAA/make_metaAA.py
#   ../make_metaAA/out/*/*.fasta