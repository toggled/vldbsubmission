#!/bin/bash
g++ -std=c++11 -g -o delmain delincore.cpp hypergraph.cpp  algorithms.cpp utils.h readhg.h

declare -a dset=("enron" "dblp")
declare -a algorithms=("Local-core-OPTIV" "deg" "clique")

del=-1 # -1 for deleting entire innermost core, if del>0 deltes del #nodes to construct h1.
writenbr=1 # 1 if want to write neighborhood dictionary as csv otherwise 0.
for dataset in "${dset[@]}"
do
    for algo in "${algorithms[@]}"
    do
        ./delmain 1 $dataset $algo $del $writenbr
        echo "------------" 
    done 
done