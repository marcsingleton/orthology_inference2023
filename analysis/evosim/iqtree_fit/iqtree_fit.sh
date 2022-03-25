# Run IQ-TREE on meta alignments and fit the rate matrix

if [ ! -d out/ ]; then
  mkdir out/
fi

../../../bin/iqtree -s ../make_meta/out/${1}.fasta -m GTR20+FO+I+R -pre out/${1} -nt $SLURM_CPUS_ON_NODE -quiet

# DEPENDENCIES
# ../make_meta/make_meta.py
#   ../make_meta/out/*.afa