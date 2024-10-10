## write a slurm script to run a job with 8 tasks on 1 node, time limit 1 hour

#!/bin/bash
#SBATCH --job-name=jobname
#SBATCH --output=jobname.out
#SBATCH --error=jobname.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=1:00:00
