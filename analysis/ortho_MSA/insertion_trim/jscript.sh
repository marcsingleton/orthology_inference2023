#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio
#SBATCH --time=12:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#
# Reporting parameters
#SBATCH --job-name=decode
#SBATCH --output=decode.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
# Link to output in scratch
if [ ! -d out ]; then
  out_dir=/global/scratch/users/singleton/orthology_inference/analysis/ortho_MSA/insertion_trim/out/
  if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}  # -p makes intermediate directory if they do not exist
  fi
  ln -s ${out_dir} out
fi

source /global/home/users/singleton/.bashrc
conda activate orthology_inference
python decode.py
