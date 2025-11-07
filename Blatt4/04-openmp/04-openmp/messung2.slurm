#!/bin/bash
#SBATCH --job-name=noel_partdiff_original
#SBATCH --open-mode=append
#SBATCH --output=ausgabe2.txt
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=west

srun make
for i in {1..10}; do
    potenz=$((2**i))
    srun ./partdiff-openmp 12 2 $potenz 2 2 16384
    echo "Lauf $i beendet"
done
echo "/n Job beendet" >> ausgabe.txt