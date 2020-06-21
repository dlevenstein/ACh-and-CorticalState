#!/bin/bash
#SBATCH -p cpu_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=120G %Increased after last fail... should also change PSS to single
echo $1
module load matlab/R2018a
matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL');addPathsDL;bz_RunAnalysis('LFPWavSpecbyDepthAnalysis','$1','basePath',true);exit;"
