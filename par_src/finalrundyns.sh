#!/bin/bash
export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=true
export OMP_PROC_BIND=true
g++ -std=c++11 -fopenmp -o optIlocalpdyn parallel_localcoreinitOptdyn.cpp    #ubuntu_cpp
# g++-11 -fopenmp -o hgmain parallel_localcore.c      #mac

declare -a dset=("enron" "bin_2" "bin_5" "congress" "contact" "dblp" "pref" "aminer")

declare -a bestThreadInit=("8" "8" "8" "8" "8" "32" "8" "32")
declare -a threads=("1" "2" "4" "8" "16" "32" "64" "128")

lb=1 # Apply Load balancing
for i in "${!dset[@]}"
    do
	data=${dset[i]}
    	bestThread=${bestThreadInit[i]}
	for t in "${threads[@]}"
	do
		./optIlocalpdyn $data $t $lb $bestThreadInit
		# echo $data.$t.$lb.$bestThread
	done
    done
