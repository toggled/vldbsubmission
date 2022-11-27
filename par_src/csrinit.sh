#!/bin/bash
export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=true
export OMP_PROC_BIND=true

g++ -std=c++11 -fopenmp -o pcsr parallel_csr.cpp    #ubuntu_cpp

#declare -a dset=("enron" "bin_2" "bin_5" "congress" "contact" "dblp" "pref" "aminer")
declare -a dset=("enron")
declare -a threads=("1" "2" "4" "8" "16" "32" "64")
#declare -a threads=("128")
for data in "${dset[@]}"
do
    for t in "${threads[@]}"
    do
	./pcsr $data $t 1
    done
done

