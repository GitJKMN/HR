#!/bin/bash
#SBATCH --job-name=original_reference_values
#SBATCH --partition=west
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=job_output.txt

# initialisation
srun ./partdiff.x 1 1 0 1 2 0
srun ./partdiff.x 1 1 0 2 2 0

# Function 1
## iterations
srun ./partdiff.x 1 1 0 1 2 1
srun ./partdiff.x 1 1 0 1 2 100
srun ./partdiff.x 1 1 100 1 2 1
srun ./partdiff.x 1 1 100 1 2 100

## precision
srun ./partdiff.x 1 1 0 1 1 1e-4
srun ./partdiff.x 1 1 0 1 1 1e-10
srun ./partdiff.x 1 1 10 1 1 1e-4
srun ./partdiff.x 1 1 10 1 1 1e-10

# Function 2
## iterations
srun ./partdiff.x 1 1 0 2 2 1
srun ./partdiff.x 1 1 0 2 2 100
srun ./partdiff.x 1 1 100 2 2 1
srun ./partdiff.x 1 1 100 2 2 100

## precision
srun ./partdiff.x 1 1 0 2 1 1e-4
srun ./partdiff.x 1 1 0 2 1 1e-10
srun ./partdiff.x 1 1 10 2 1 1e-4
srun ./partdiff.x 1 1 10 2 1 1e-10

# initialisation
srun ./partdiff.x 1 2 0 1 2 0
srun ./partdiff.x 1 2 0 2 2 0

# Function 1
## iterations
srun ./partdiff.x 1 2 0 1 2 1
srun ./partdiff.x 1 2 0 1 2 100
srun ./partdiff.x 1 2 100 1 2 1
srun ./partdiff.x 1 2 100 1 2 100

## precision
srun ./partdiff.x 1 2 0 1 1 1e-4
srun ./partdiff.x 1 2 0 1 1 1e-10
srun ./partdiff.x 1 2 10 1 1 1e-4
srun ./partdiff.x 1 2 10 1 1 1e-10

# Function 2
## iterations
srun ./partdiff.x 1 2 0 2 2 1
srun ./partdiff.x 1 2 0 2 2 100
srun ./partdiff.x 1 2 100 2 2 1
srun ./partdiff.x 1 2 100 2 2 100

## precision
srun ./partdiff.x 1 2 0 2 1 1e-4
srun ./partdiff.x 1 2 0 2 1 1e-10
srun ./partdiff.x 1 2 10 2 1 1e-4
srun ./partdiff.x 1 2 10 2 1 1e-10