# Run PhyML on meta alignments using GTR model

if [ ! -d out/${1}/ ]; then
  mkdir -p out/${1}/
fi

../../../bin/phyml -i ../make_meta/out/${1}/${2}.phy --model GTR --pinv ${3} -o tlr --search SPR --rand_start --n_rand_starts 3 > out/${1}/${2}.out
mv ../make_meta/out/${1}/${2}.phy_phyml_rand_trees.txt ../make_meta/out/${1}/${2}.phy_phyml_stats.txt ../make_meta/out/${1}/${2}.phy_phyml_tree.txt out/${1}/

# DEPENDENCIES
# ../make_meta/make_meta.py
#   ../make_meta/out/*/*.phy