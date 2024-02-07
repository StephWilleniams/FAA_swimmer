#!/bin/bash

#SBATCH --nodes=1  #asked for 1 node
#SBATCH --ntasks=10 #asked for 1 cores
# #SBATCH --partition medium  #this job will submit to medium partition
#SBATCH --array=1-5
#SBATCH --mem=1G  #this job is asked for 1G of total memory, use 0 if you want to use entire node memory
#SBATCH --time=0-00:15:00 # 15 minutes
#SBATCH --output=test_%A_%a.qlog  #the output information will put into test_$SLURM_ARRAY_JOB_ID_$SLURM_ARRAY_TASK_ID.qlog file
#SBATCH --job-name=test1  #the job name
#SBATCH --export=ALL