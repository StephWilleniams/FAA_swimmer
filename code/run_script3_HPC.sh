#!/bin/bash
#SBATCH --nodes=1  #asked for 1 node
#SBATCH --ntasks=40 #asked for 1 cores
#SBATCH --array=1-12
#SBATCH --mem=1G  #this job is asked for 1G of total memory, use 0 if you want to use entire node memory
#SBATCH --time=0-06:00:00 # 15 minutes
#SBATCH --output=test_%A_%a.qlog  #the output information will put into test_$SLURM_ARRAY_JOB_ID_$SLURM_ARRAY_TASK_ID.qlog file
#SBATCH --job-name=hiQueueViewer  #the job name
#SBATCH --export=ALL

Nruns=40

# Iterate through numbers 1 to 5
for i in {1..$Nrun} # Does it do endpoint??
do
  mkdir outputs/${SLURM_ARRAY_TASK_ID}_${i}_outputs/
  /a.out 20 -1 $SLURM_ARRAY_TASK_ID $i &
done

mkdir outputs/${SLURM_ARRAY_TASK_ID}_${Nruns}_outputs/
./a.out $SLURM_ARRAY_TASK_ID $Nruns 

wait
