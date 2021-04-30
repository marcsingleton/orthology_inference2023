#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio2
#SBATCH --time=18:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#
# Reporting parameters
#SBATCH --job-name=feature_calc_shuffle
#SBATCH --output=feature_calc_shuffle.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
source /global/home/users/singleton/.bashrc
conda activate IDREvoDevo
python feature_calc_shuffle.py
