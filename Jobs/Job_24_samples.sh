#!/bin/bash
#SBATCH --job-name=FRASER_24_samples
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-user="j.m.van_der_molen@lumc.nl"
#SBATCH --mail-type="ALL"
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

echo "loading modules"
module purge
module load tools/miniconda/python3.12/24.9.2

echo "activating conda environment"
eval "$(conda shell.bash hook)"
conda activate /exports/sascstudent/jmvandermolen1/conda_envs/drop_env

cd /exports/sascstudent/jmvandermolen1/Developer/drop_demo_fraser2

echo "Running FRASER on 24 samples"

snakemake aberrantSplicing \
    --cores "$SLURM_CPUS_PER_TASK" \
    --rerun-triggers mtime \
    -p

echo "FRASER completed!"