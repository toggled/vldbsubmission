#!/bin/bash
g++ -std=c++11 -g -o delmain delincore.cpp hypergraph.cpp  algorithms.cpp utils.h readhg.h
mkdir -p "../python_src/sirdata_naheed_vai"

declare -a dset=("enron")
# declare -a dset=("pref" "aminer") 
declare -a algorithms=("Local-core-OPTIV" "deg" "clique")

# del=-1 # -1 for deleting entire innermost core, if del>0 deltes del #nodes to construct h1.
# writenbr=1 # 1 if want to write neighborhood dictionary as csv otherwise 0.
# for dataset in "${dset[@]}"
# do
#     for algo in "${algorithms[@]}"
#     do
#         ./delmain 1 $dataset $algo $del $writenbr
#         echo "------------" 
#     done 
# done

del=10
writenbr=1 # 1 if want to write neighborhood dictionary as csv otherwise 0.
for dataset in "${dset[@]}"
do
    for algo in "${algorithms[@]}"
    do
        ./delmain 1 $dataset $algo $del $writenbr
        echo "------------" 
    done 
done
del=5
writenbr=1 # 1 if want to write neighborhood dictionary as csv otherwise 0.
for dataset in "${dset[@]}"
do
    for algo in "${algorithms[@]}"
    do
        ./delmain 1 $dataset $algo $del $writenbr
        echo "------------" 
    done 
done
