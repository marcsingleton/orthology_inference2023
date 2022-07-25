# Run IQ-TREE on meta alignments using GTR model

if [ ! -d out/${1}/ ]; then
  mkdir -p out/${1}/
fi

../../../bin/iqtree -s ../make_meta_NT/out/${1}/${2}.afa -m GTR+FO${3}+G -pre out/${1}/${2} -quiet

# DEPENDENCIES
# ../make_meta_NT/make_meta_NT.py
#   ../make_meta_NT/out/*/*.afa