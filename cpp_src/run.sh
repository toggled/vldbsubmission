#!/bin/bash

# g++-11 -Wall -g  -o main main.cpp hypergraph.cpp  algorithms.cpp readhg.h utils.h #  run on clang
g++ -std=c++11 -Wall -g -o main main.cpp hypergraph.cpp  algorithms.cpp utils.h readhg.h #  run on gnu c++ compiler


declare -a dset=("enron")
#  "bin_1" "bin_2" "bin_4" "bin_5"  "contact" "congress" "dblp", "aminer")
declare -a algorithms=("Peel" "E-Peel" "Local-core" "Local-core-OPTI" "Local-core-OPTII" "Local-core-OPTIII" "Local-core-OPTIV")

it=5 # #Iterations to run each algorithm on each dataset.
log=0 # Activate logging to output core-numbers & iteration h-index statistics
for dataset in "${dset[@]}"
do
    for algo in "${algorithms[@]}"
    do
        ./main 1 $dataset $algo $it $log
        # echo "------------" 
    done 
done