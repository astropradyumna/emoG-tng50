#!/bin/bash -l
#SBATCH -J FoFRad
#SBATCH -p saleslab
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=7
#SBATCH --mem=300gb
###SBATCH --mem-per-cpu=G
#SBATCH --time=48:00:00
#SBATCH -o output_log/epir1.out
#SBATCH -e output_log/epir1.err
#SBATCH --mail-user=psadh003@ucr.edu
#SBATCH --mail-type=ALL

# Load needed modules
# You could also load frequently used modules from within your ~/.bashrc
module load slurm # Should already be loaded
module load openmpi # Should already be loaded
#module load hdf5

cd /rhome/psadh003/bigdata/tng50/data_extractors

python3 extract_particles_inside_rad.py 2