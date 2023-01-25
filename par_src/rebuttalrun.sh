#!/bin/bash
export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=true
export OMP_PROC_BIND=true

g++ -std=c++11 -fopenmp -o optIlocalp parallel_localcoreinitOpt2.cpp
# g++-11 -fopenmp -o hgmain parallel_localcore.c      #mac
#declare -a dset=("enron" "bin_2" "bin_5" "contact")
declare -a dset=("pref3U_50mil_21_1_500")
#declare -a dset=("congress" "dblp" "pref")

declare -a threads=("32" "64" "128")

lb=1 # Apply Load balancing
bestThread=64
for data in "${dset[@]}"
    do
    for t in "${threads[@]}"
        do
		./optIlocalp $data $t $lb $bestThread
        done
    done 
