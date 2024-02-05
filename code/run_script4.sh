#!/bin/bash

# # Tidy directories
rm outputs/*/*
rmdir outputs/*
mkdir outputs

# # module load parallel
g++ -O3 main.cpp # Compile main, pre-prepare this
Nruns=1

# # Iterate through numbers 1 to Nruns
for (( i=1; i<=$Nruns; i++ )) # Does it do endpoint??                                                                                                                            
do

  mkdir "outputs/1_${i}_outputs/"
  ./a.out 20 -1 1 $i 

done
