#!/bin/bash
#SBATCH --account=nn11072k
#SBATCH --job-name=Rumigen-paper-1-1.005.5.pg
#SBATCH --partition=normal
#SBATCH --time=0-00:05:00
#SBATCH --ntasks=20
#SBATCH --mem=64G

# it is good to have the following lines in any bash script
set -o errexit  # make bash exit on any error
set -o nounset  # treat unset variables as errors

#cost # to show how many hours left for this project
#dusage -p nn11072k # to show how much storage left for this project
./