#!/bin/bash
#
#SBATCH --job-name=juliamicro
#SBATCH --output=./test_script_outputs/output%A_%a.txt
#
#SBATCH --ntasks=2
#SBATCH --array=0-15
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --exclude=node128


#srun julia microGenerator.jl IEEE-RTS96.csv $SLURM_ARRAY_TASK_ID
echo $SLURM_JOBID
echo $SLURM_ARRAY_TASK_ID
