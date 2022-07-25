# Run IQ-TREE on meta alignments using LG model

if [ ! -d out/${1}/ ]; then
  mkdir -p out/${1}/
fi

../../../bin/iqtree -s ../make_meta_AA/out/${1}/${2}.afa -m LG+FO${3}+G -pre out/${1}/${2} -quiet

# DEPENDENCIES
# ../make_meta_AA/make_meta_AA.py
#   ../make_meta_AA/out/*/*.afa