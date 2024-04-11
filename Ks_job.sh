#!/bin/bash
#SBATCH -t 6-00:00
#SBATCH --job-name=Kstar
#SBATCH --partition=htc
#SBATCH --cluster=smp
#SBATCH --mail-user=zakaria@pitt.edu
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100g
#SBATCH --cpus-per-task=28

module load gcc/12.2.0
module load python/ondemand-jupyter-python3.8
module load intel/2017.1.132

pip3 install pandas
pip3 install numpy
pip3 install statsmodels

# first arg is path to rds and second is path where you want to sve results
python3 ./SIV-research-Tools/exce.py 28 ./SIV_SGS_files_from_Ambrose-Lin_study/V1_V5/gp120_annotated_dfs/ ./SIV_SGS_files_from_Ambrose-Lin_study/V1_V5/gp120_annotated_dfs/

