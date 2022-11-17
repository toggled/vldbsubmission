#!/bin/bash

# g++-11 -Wall -g  -o kdmain kdmain.cpp hypergraph.cpp  algorithms.cpp readhg.h utils.h #  run on clang
g++ -std=c++11 -Wall -g -o kdmain kdmain.cpp hypergraph.cpp  algorithms.cpp utils.h readhg.h #  run on gnu c++ compiler


declare -a dset=("congress" "enron" "dblp" "pref" "aminer")
declare -a algorithms=("kdcore")

it=1 # #Iterations to run each algorithm on each dataset.
log=0 # logging is not implemented
for dataset in "${dset[@]}"
do
    for algo in "${algorithms[@]}"
    do
        ./kdmain 1 $dataset $algo $it
        echo "------------" 
    done 
done

# output order (results.csv): algo,dataset,execution time,init_time,num_threads,total iteration(p),total iteration(s)
# output order (kdcore.csv): node id, primary core, secondary core.
