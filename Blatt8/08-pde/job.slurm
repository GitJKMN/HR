#!/bin/bash
#SBATCH --job-name=partdiff_mpi
#SBATCH --output=partdiff_%j.out
#SBATCH --error=STDOUT
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=west
#SBATCH --time=00:30:00

echo "=== RUN 1: 1 Prozess (sequenziell) ==="
mpirun -np 1 ./partdiff.x 1 2 512 2 2 800

echo "=== RUN 2: 5 Prozesse (parallel) ==="
mpirun -np 5 ./partdiff.x 1 2 512 2 2 800

echo "=== RUN 3: 12 Prozesse (parallel) ==="
mpirun -np 12 ./partdiff.x 1 2 512 2 2 800
