#!/bin/bash

# g++-11 -Wall -g  -o main main.cpp hypergraph.cpp  algorithms.cpp readhg.h utils.h #  run on clang
g++ -std=c++11 -Wall -g -o main2 main.cpp hypergraph.cpp  algorithms.cpp utils.h readhg.h #  run on gnu c++ compiler


#declare -a dset=("bin_2" "bin_5"  "contact" "congress" "enron" "dblp" "pref" "aminer")
declare -a dset=("aminer")
declare -a algorithms=("Peel")

it=2 # #Iterations to run each algorithm on each dataset.
log=1 # Activate logging to output core-numbers & iteration h-index statistics
for algo in "${algorithms[@]}"
do
    for dataset in "${dset[@]}"
    do
        ./main2 1 $dataset $algo $it $log
        echo "------------" 
    done 
done
