
#!/bin/bash

g++ -O3 main.cpp

rm 1_outputs/*
rmdir 1_outputs/
mkdir 1_outputs/
./a.out 0.2 1 1

rm 2_outputs/*
rmdir 2_outputs/
mkdir 2_outputs/
./a.out 5.6 1 2
