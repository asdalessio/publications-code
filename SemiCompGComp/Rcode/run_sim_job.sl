#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=40g
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 5-00:00
#SBATCH --array 1-2
#SBATCH --output=/dev/null --error=/dev/null
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=asoren@email.unc.edu

# module purge
module add r/4.4.0

Rscript 030326_multistate_g_comp_simulation_publish.R $SLURM_ARRAY_TASK_ID
