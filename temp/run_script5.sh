#!/bin/bash                                                                                                                                                                     
#SBATCH --nodes=1  #asked for 1 node                                                                                                                                            
#SBATCH --ntasks=40 #asked for 1 cores                                                                                                                                          
# #SBATCH --partition medium  #this job will submit to medium partition                                                                                                         
# #SBATCH --array=1-12                                                                                                                                                             
#SBATCH --mem=1G  #this job is asked for 1G of total memory, use 0 if you want to use entire node memory                                                                        
#SBATCH --time=0-00:15:00 # 15 minutes                                                                                                                                          
#SBATCH --output=test_%A_%a.qlog  #the output information will put into test_$SLURM_ARRAY_JOB_ID_$SLURM_ARRAY_TASK_ID.qlog file                                                 
#SBATCH --job-name=test1  #the job name                                                                                                                                         
#SBATCH --export=ALL      

# # Tidy directories

# # module load parallel
# g++ -O3 main.cpp # Compile main
Nruns=10
runind=0

# # Iterate through numbers 1 to Nruns
for (( i=1; i<=$Nruns; i++ )) # Does it do endpoint??                                                                                                                            
do
  runind=$(($runind+1))
  mkdir "outputs/1_1_${runind}_outputs/" 
  ./a.out 20 -1 1 1 $runInd &
done

for (( i=1; i<=$Nruns; i++ )) # Does it do endpoint??                                                                                                                            
do
  runind=$(($runind+1))
  mkdir "outputs/12_1_${runind}_outputs/"  
  ./a.out 20 -1 12 1 &
done

for (( i=1; i<=$Nruns; i++ )) # Does it do endpoint??                                                                                                                            
do
  runind=$(($runind+1))
  mkdir "outputs/1_40_${runind}_outputs/"
  ./a.out 20 -1 1 40 &
done

for (( i=1; i<$Nruns; i++ )) # Does it do endpoint??                                                                                                                            
do
  runind=$(($runind+1))
  mkdir "outputs/12_40_${runind}_outputs/"
  ./a.out 20 -1 12 40 &
done

runind=$(($runind+1))
mkdir "outputs/12_40_${runind}_outputs/"
./a.out 20 -1 $SLURM_ARRAY_TASK_ID 4 

wait
