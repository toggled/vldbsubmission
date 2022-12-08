#!/bin/bash
g++ -std=c++11 -Wall -g -o bpmain main.cpp hypergraph.cpp  algorithms.cpp utils.h readhg.h 
declare -a dset=("bin_2" "bin_5"  "contact" "congress" "enron" "dblp" "pref" "aminer")
for dataset in "${dset[@]}"
do
	./bpmain 1 $dataset bipartite 1 1

done

