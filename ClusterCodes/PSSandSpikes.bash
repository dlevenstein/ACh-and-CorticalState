#!/bin/bash
#SBATCH -p cpu_short
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=32G
echo $1
module load matlab/R2018a
matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/ClusterCodes');addPaths;cd('/gpfs/data/buzsakilab/DL');addPathsDL;bz_RunAnalysis('PSSandSpikesAnalysis','$1','basePath',true);exit;"
