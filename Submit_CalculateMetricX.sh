#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 12
#SBATCH --mem 8000

#SBATCH --job-name=CalculateMetric.job

#SBATCH --output=/home/swenson/projects/PEX_feature_metrics/CalculateMetric.out

#SBATCH --error=/home/swenson/projects/PEX_feature_metrics/CalculateMetric.err

#SBATCH -p med

#SBATCH --time=6:00:00

source activate /home/swenson/.bashrc

conda activate heatwave

python3 /home/swenson/projects/PEX_feature_metrics/CalculateVortical_Metric.py > CalculateMetric.log