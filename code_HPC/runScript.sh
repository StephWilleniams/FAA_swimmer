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

Nruns=10 # Equal to the ntasks value                                                                                                                                            

# Iterate through numbers 1 to 5                                                                                                                                                
for (( i=1; i<$Nruns; i++ )) # Does it do endpoint??                                                                                                                            
do

    rm "${SLURM_ARRAY_TASK_ID}_${i}_outputs/*"
    rmdir "${SLURM_ARRAY_TASK_ID}_${i}_outputs/"
    mkdir "${SLURM_ARRAY_TASK_ID}_${i}_outputs/"
    # echo  "$SLURM_ARRAY_TASK_ID $i" # Test outputs                                                                                                                            
    /a.out 20 -1 $SLURM_ARRAY_TASK_ID $i &

done

rm "${SLURM_ARRAY_TASK_ID}_${Nruns}_outputs/*"
rmdir "${SLURM_ARRAY_TASK_ID}_${Nruns}_outputs/"
mkdir "${SLURM_ARRAY_TASK_ID}_${Nruns}_outputs/"
# echo  "$SLURM_ARRAY_TASK_ID $Nruns" # Test outputs                                                                                                                            
./a.out 20 -1 $SLURM_ARRAY_TASK_ID $Nruns

wait