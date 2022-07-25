# Run IQ-TREE on meta alignments using GTR2 model

if [ ! -d out/${1}/ ]; then
  mkdir -p out/${1}/
fi

../../../bin/iqtree -s ../make_meta_INDEL/out/${1}/${2}.afa -m GTR2+FO${3} -pre out/${1}/${2} -quiet

# DEPENDENCIES
# ../make_meta_INDEL/make_meta_INDEL.py
#   ../make_meta_INDEL/out/*/*.afa