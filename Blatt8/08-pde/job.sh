#!/bin/bash
#SBATCH --job-name=original_reference_values
#SBATCH --partition=west
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=job_output2.txt
#SBATCH --error=job_output2.txt

# Slurm gibt die Anzahl Tasks
NTASKS=$SLURM_NTASKS

# ----------------------------
# Alle Tests nacheinander, jeder Test parallel auf 2 Prozessen
# ----------------------------
TESTS=(
  "1 1 0 1 2 0"
  "1 1 0 2 2 0"
  "1 1 0 1 2 1"
  "1 1 0 1 2 100"
  "1 1 100 1 2 1"
  "1 1 100 1 2 100"
  "1 1 0 1 1 1e-4"
  "1 1 0 1 1 1e-10"
  "1 1 10 1 1 1e-4"
  "1 1 10 1 1 1e-10"
  "1 1 0 2 2 1"
  "1 1 0 2 2 100"
  "1 1 100 2 2 1"
  "1 1 100 2 2 100"
  "1 1 0 2 1 1e-4"
  "1 1 0 2 1 1e-10"
  "1 1 10 2 1 1e-4"
  "1 1 10 2 1 1e-10"
  "1 2 0 1 2 0"
  "1 2 0 2 2 0"
  "1 2 0 1 2 1"
  "1 2 0 1 2 100"
  "1 2 100 1 2 1"
  "1 2 100 1 2 100"
  "1 2 0 1 1 1e-4"
  "1 2 0 1 1 1e-10"
  "1 2 10 1 1 1e-4"
  "1 2 10 1 1 1e-10"
  "1 2 0 2 2 1"
  "1 2 0 2 2 100"
  "1 2 100 2 2 1"
  "1 2 100 2 2 100"
  "1 2 0 2 1 1e-4"
  "1 2 0 2 1 1e-10"
  "1 2 10 2 1 1e-4"
  "1 2 10 2 1 1e-10"
)

for test in "${TESTS[@]}"; do
  echo "Running: ./partdiff.x $test"
  mpiexec -n $NTASKS ./partdiff.x $test
done
