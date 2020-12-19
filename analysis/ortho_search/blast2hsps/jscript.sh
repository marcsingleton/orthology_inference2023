#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio2
#SBATCH --time=8:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#
# Reporting parameters
#SBATCH --job-name=blast2hsps
#SBATCH --output=blast2hsps.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
if [ ! -d out ]; then
  out_dir=/global/scratch/singleton/IDREvoDevo/analysis/ortho_cluster2/blast2hsps/out
  if [ ! -d $out_dir ]; then
    mkdir -p $out_dir  # -p makes intermediate directory if they do not exist
  fi
  ln -s $out_dir out
fi

module load python
python blast2hsps.py
