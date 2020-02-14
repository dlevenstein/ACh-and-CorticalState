#!/bin/bash
#SBATCH -p cpu_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=20G
echo $1
module load matlab/R2018a
matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL');addPathsDL;bz_RunAnalysis('BehaviorAnalysis2','$1','basePath',true);exit;"
