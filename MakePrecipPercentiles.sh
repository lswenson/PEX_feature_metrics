#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 32000

#SBATCH --job-name=MakePrecipPercentiles.job

#SBATCH --output=/home/swenson/projects/PEX_feature_metrics/MakePrecipPercentiles.out

#SBATCH --error=/home/swenson/projects/PEX_feature_metrics/MakePrecipPercentiles.err

#SBATCH -p med

#SBATCH --time=24:00:00

source activate /home/swenson/.bashrc

conda activate heatwave

python3 /home/swenson/projects/PEX_feature_metrics/MakePrecipPercentiles.py > MakePrecipPercentiles.log