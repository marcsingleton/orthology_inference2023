#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio2
#SBATCH --time=48:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#
# Reporting parameters
#SBATCH --job-name=iqtree_fit
#SBATCH --output=iqtree_fit.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
source /global/home/users/singleton/.bashrc
conda activate orthology_inference
module load gnu-parallel

for file in $(ls ../make_meta/out/); do
  if [[ $file == *.afa ]]; then
    bash iqtree_fit.sh $(basename $file .afa)
  fi
done
