#!/bin/bash
#SBATCH --job-name=DROP_smoketest_2
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module purge
module load tools/miniconda/python3.9/4.12.0
eval "$(conda shell.bash hook)"
conda activate /exports/sascstudent/jmvandermolen1/conda_envs/drop_env

cd /exports/sascstudent/jmvandermolen1/Developer/drop_demo_fraser2

echo "=== Dry-run check ==="
snakemake aberrantSplicing --cores 1 -n -r