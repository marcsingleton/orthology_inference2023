#!/bin/bash
# Run job script, substituting species name and path as arguments

# Set current directory and link to output in scratch
cd "$(dirname "$0")"
if [ ! -d out ]; then
  out_dir=/global/scratch/users/singleton/IDREvoDevo/analysis/ortho_search/blast_search/out/
  if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}  # -p makes intermediate directory if they do not exist
  fi
  ln -s ${out_dir} out
fi

# Queue up commands
cat ../config/genomes.tsv | while read spid txid source prot_path tcds_path
do
  if [[ ${spid} != \#* ]]  # Double brackets is expanded syntax for tests
  then
    sbatch << _EOF_
#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio2
#SBATCH --time=48:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#
# Reporting parameters
#SBATCH --job-name=blast_search_${spid}
#SBATCH --output=out/blast_search_${spid}.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
source /global/home/users/singleton/.bashrc
conda activate IDREvoDevo
python blast_search.py ${spid} ../remove_duplicates/out/${spid}.fa ../../../bin/ncbi-blast-2.13.0+/bin/blastp
_EOF_
  fi
done