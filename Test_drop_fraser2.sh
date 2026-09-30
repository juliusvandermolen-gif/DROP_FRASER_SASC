#!/bin/bash
#SBATCH --job-name=drop_splicing
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module purge
module load tools/miniconda/python3.9/4.12.0
eval "$(conda shell.bash hook)"
conda activate /exports/sascstudent/jmvandermolen1/conda_envs/drop_env
module list

snakemake aberrantSplicing \
    --cores $SLURM_CPUS_PER_TASK \
    --rerun-triggers mtime \