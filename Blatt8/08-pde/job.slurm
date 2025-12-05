#!/bin/bash
#SBATCH --job-name=partdiff3
#SBATCH --output=partdiff_%j.out
#SBATCH --error=STDOUT
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=west

echo "=== RUN 1 ==="
srun --ntasks=1 ./partdiff.x 1 2 512 2 2 800

echo "=== RUN 2 ==="
srun --ntasks=5 ./partdiff.x 1 2 512 2 2 800

echo "=== RUN 3 ==="
srun --ntasks=12 ./partdiff.x 1 2 512 2 2 800
