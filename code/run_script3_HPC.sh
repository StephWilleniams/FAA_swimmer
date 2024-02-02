#!/bin/bash  
#SBATCH --nodes=1  # asked for 1 node
#SBATCH --ntasks=20 # asked for 20 cores
#SBATCH --partition test  # this job will submit to test partition
#SBATCH --mem=96G  #this job is asked for 96G of total memory, use 0 if you want to use entire node memory
# #SBATCH --gres=gpu:X # uncomment this line if you need GPU access, replace X with number of GPU you need
# #SBATCH -w <selected_node> #uncomment this line if you want to select specific available node to run 
#SBATCH --time=0-00:15:00 # 15 minutes  
#SBATCH --output=test1.qlog  #the output information will put into test1.qlog file
#SBATCH --job-name=test1  #the job name
#SBATCH --export=ALL

module load parallel

g++ -O3 main.cpp # Compile main

kmax = 10 # max kick frequency to calculate

srun="srun -n1 -N1 --exclusive"

parallel="parallel -N 1 --delay .2 -j $SLURM_NTASKS --joblog parallel_joblog --resume"

$parallel "$srun ./runtask arg1:{1}" ::: {1..40}
do

    kf = 100*$i*$kmax/40
    ks = 100*1
    echo "$kf"
    echo "$ks"
    ./a.out 20 -1 $kf $ks $i # Arguments: (100*dR,polarity,100*kickFreq,100*kickStr,fileLabel)

done
