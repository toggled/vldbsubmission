#!/bin/bash

g++ -std=c++11 -Wall -g -o densmain densest_subhypergraph.cpp hypergraph.cpp readhg.h #  run on gnu c++ compiler

# declare -a dset=("contact")
declare -a dset=("klay" "protein" "bin_2" "bin_5"  "contact" "congress" "enron" "dblp" "pref" "aminer")
declare -a algorithms=("nbr" "deg")

for dataset in "${dset[@]}"
do
    for algo in "${algorithms[@]}"
    do
        ./densmain $dataset $algo 
        echo "------------" 
    done 
done