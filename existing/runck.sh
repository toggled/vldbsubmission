g++ -std=c++11 -g -o ckmain clique_graph.cpp hypergraph.cpp  algorithms.cpp readhg.h
#./ckmain 1 aminer clique
#g++ -std=c++11 -g -o bpmain bipartiteGraph.cpp hypergraph.cpp  algorithms.cpp readhg.h
#./bpmain 1 aminer bipartite
declare -a dset=("bin_2" "bin_5"  "contact" "congress" "enron" "dblp" "pref" "aminer")
for dataset in "${dset[@]}"
do
	./ckmain 1 $dataset clique

done

