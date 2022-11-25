#g++ -std=c++11 -g -o ckmain clique_graph.cpp hypergraph.cpp  algorithms.cpp readhg.h
#./ckmain 1 aminer clique
g++ -std=c++11 -g -o bpmain bipartiteGraph.cpp hypergraph.cpp  algorithms.cpp readhg.h
./bpmain 1 aminer bipartite
