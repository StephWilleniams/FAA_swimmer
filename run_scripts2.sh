
#!/bin/bash

g++ -O3 main.cpp # Compile main

# Arguments: (100*dR,polarity,100*kickFreq,100*kickStr,fileLabel)

# Puller, no kicks
# rm 1_outputs/*
# rmdir 1_outputs/
# mkdir 1_outputs/
# ./a.out 20 -1 0 0 1

# Pusher, no kicks
# rm 2_outputs/*
# rmdir 2_outputs/
# mkdir 2_outputs/
# ./a.out 560 1 0 0 2

# Puller, strong/likely kicks
# rm 3_outputs/*
# rmdir 3_outputs/
# mkdir 3_outputs/
# ./a.out 20 -1 100 100 3

# Pusher, strong/likely kicks
rm 4_outputs/*
rmdir 4_outputs/
mkdir 4_outputs/
./a.out 560 1 100 100 4
