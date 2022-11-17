#!/bin/bash
export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=true
export OMP_PROC_BIND=true

#g++ -std=c++11 -fopenmp -o hgmain parallel_localcore.cpp    #ubuntu_cpp
# g++-11 -fopenmp -o hgmain parallel_localcore.c      #mac

declare -a dset=("enron" "bin_2" "bin_5" "congress" "contact" "dblp" "pref" "aminer")

declare -a threads=("1" "2" "4" "8" "16" "32" "64")
#declare -a threads=("128")
lb=1 # Apply Load balancing
for data in "${dset[@]}"
    do
    for t in "${threads[@]}"
        do
            ./hgmain $data $t $lb
        done
    done 

lb=0 # No load-balancing
for data in "${dset[@]}"
    do
    # declare -a threads=("1" "2" "4" "8")
    for t in "${threads[@]}"
        do
            ./hgmain $data $t $lb
        done
    done 
