#!/bin/bash
#SBATCH --job-name=noel_partdiff_original
#SBATCH --open-mode=append
#SBATCH --output=ausgabe.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=west

srun make
for i in {1..12}; do
    srun ./partdiff-openmp $i 2 512 2 2 1024
done
echo "Job beendet" >> ausgabe.txt
