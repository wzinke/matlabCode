#!/bin/bash
#SBATCH --mail-user=kaleb.a.lowe@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=500M
#SBATCH --output=matlab_job_slurm.out

setpkgs -a matlab_r2016a

matlab -nodisplay -nosplash < ./testNewDir/testSlurm.m
