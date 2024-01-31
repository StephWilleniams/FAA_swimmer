#!/bin/bash

g++ -O3 main.cpp # Compile main

kmax = 10 # max kick frequency to calculate

for i in {1..16}
do

    kf = 100*$i*$kmax/16
    ks = 100*1
    echo "$kf"
    echo "$ks"
    ./a.out 20 -1 $kf $ks $i # Arguments: (100*dR,polarity,100*kickFreq,100*kickStr,fileLabel)

done
