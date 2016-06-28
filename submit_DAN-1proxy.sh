#!/bin/bash
#
#SBATCH --job-name=juliaDAN-1proxy
#SBATCH --output=output.txt
#
#SBATCH --array=0-1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=1000


srun julia DAN-1proxy.jl IEEE-RTS96.csv $SLURM_ARRAY_TASK_ID

