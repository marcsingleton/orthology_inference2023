# Run IQ-TREE on meta alignments using WAG model

if [ ! -d out/${1}/ ]; then
  mkdir -p out/${1}/
fi

../../../bin/iqtree -s ../make_AAmeta/out/${1}/${2}.fasta -m WAG+FO${3}+G -pre out/${1}/${2} -quiet

# DEPENDENCIES
# ../make_meta/make_AAmeta.py
#   ../make_AAmeta/out/*/*.phy