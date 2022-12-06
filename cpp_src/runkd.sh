#!/bin/bash

# g++-11 -Wall -g  -o main main.cpp hypergraph.cpp  algorithms.cpp readhg.h utils.h #  run on clang
g++ -std=c++11 -o kdmain main.cpp hypergraph.cpp  algorithms.cpp utils.h readhg.h #  run on gnu c++ compiler


declare -a dset=("congress" "enron" "dblp" "pref" "aminer")
declare -a algorithms=("kdcore")

it=2 # #Iterations to run each algorithm on each dataset.
log=0 # Activate logging to output core-numbers & iteration h-index statistics
for dataset in "${dset[@]}"
do
    for algo in "${algorithms[@]}"
    do
        ./kdmain 1 $dataset $algo $it $log
        echo "------------" 
    done 
done

# declare -a dset=("klay" "protein")
# declare -a algorithms=("kdcore")

# it=1 # #Iterations to run each algorithm on each dataset.
# log=0 # Activate logging to output core-numbers & iteration h-index statistics
# for dataset in "${dset[@]}"
# do
#     for algo in "${algorithms[@]}"
#     do
#         ./main 1 $dataset $algo $it $log
#         echo "------------" 
#     done 
# done
