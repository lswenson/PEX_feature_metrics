#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 12
#SBATCH --mem 16000

#SBATCH --job-name=MakeWetDays.job

#SBATCH --output=/home/swenson/projects/PEX_feature_metrics/MakeWetDays.out

#SBATCH --error=/home/swenson/projects/PEX_feature_metrics/MakeWetDays.err

#SBATCH -p med

#SBATCH --time=18:00:00

source activate /home/swenson/.bashrc

conda activate heatwave

python3 /home/swenson/projects/PEX_feature_metrics/MakeWetDays.py > MakeWetDays.log