#!/bin/bash
#SBATCH --job-name=plot_pdf_FRASER
#SBATCH --partition="short"
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-user="j.m.van_der_molen@lumc.nl"
#SBATCH --mail-type="ALL"
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module purge
module load tools/miniconda/python3.9/4.12.0
eval "$(conda shell.bash hook)"
conda activate /exports/sascstudent/jmvandermolen1/conda_envs/drop_env
module list
 
Rscript Plot_R.R