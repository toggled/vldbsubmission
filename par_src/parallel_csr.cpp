#define _GNU_SOURCE
#include <omp.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <string> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <queue>
#include <functional> 
#include <cassert>
#include <unordered_map>
typedef  std::vector<std::string> strvec;
typedef  std::vector<size_t> intvec;
typedef  std::set<size_t> intset;
typedef std::map<std::string, std::string> strstrMap;
// typedef std::map<size_t, intvec> intintVMap;
typedef std::map<size_t, intset> intintSMap;
typedef std::vector< intvec > intintVec;
typedef std::pair<size_t, size_t> pi;
typedef std::unordered_map<size_t,size_t> intintMap;

#ifdef DEBUG
#  define D(x) x
#else
#  define D(x) 
#endif

#ifdef OUTPUT
#  define O(x) x
#else
#  define O(x) 
#endif
	
strstrMap dataset_to_filename = {
 			{"enron" , "../data/datasets/real/Enron.hyp"},
            {"congress" , "../data/datasets/real/congress-bills.hyp"},
            {"contact" , "../data/datasets/real/contact-primary-school.hyp"},
            {"dblp", "../data/datasets/real/DBLP.hyp"},
            {"aminer","../data/datasets/real/aminer.hyp"},
            {"bin_2" , "../data/datasets/synthetic/binomial_5_500_4_0.200000_sample_2_iter_1.txt"},
            {"bin_5" , "../data/datasets/synthetic/binomial_5_500_3_0.200000_sample_5_iter_1.txt"},
            {"pref", "../data/datasets/synthetic/pref_1000000_3_1.hyp"}
        };

typedef struct edge {
	/* Edge container.
	 * Also doubles as a message.
	 */
	size_t id;
	size_t est;
	struct edge *next;
} edge;

typedef struct graph_node {
	/* A vertex in a graph.
	 */
	size_t id;
	size_t valid;
	size_t kcore;
	size_t active;
	size_t num_nbrs;
	std::vector<int> h_count;

} graph_node;

void re_order_loadbalance( size_t n, size_t k, intintVec& B, intvec& prefixsum_B, intvec& node_index_index, intintVec& hyperedges,intintMap& node_index, size_t nbrs_N[]){
	std::priority_queue<pi, std::vector<pi>, std::greater<pi> > pq;
	for (size_t i = 0; i<k; i++){
		std::vector<size_t> partitionA;
		B.push_back(partitionA);
		pq.push(std::make_pair(0, i));
	}
	for (size_t i = 0; i<n; i+=k){
		for(size_t j=i; j< std::min(n,i+k); j++){
			std::pair<size_t,size_t> max_pair = pq.top();
			B[max_pair.second].push_back(j);
			max_pair.first += (nbrs_N[j+1]-nbrs_N[j]);
			pq.pop();
			pq.push(max_pair);
		}
	}
}
void no_loadbalance(size_t n, size_t k, intintVec& B, intvec& prefixsum_B, intvec& node_index_index, intintVec& hyperedges,intintMap node_index, size_t nbrs_N[]){
 	size_t  i,j,kk;
 	std::cout<<"simple loadbalance\n";
 	for (size_t i = 0; i<k; i++){
 		std::vector<size_t> partitionA;
 		B.push_back(partitionA);
 	}
 	for (i = 0,j=0; i<n; i+=n/k,j++){
 		if (j<k){
 			for(kk = i; kk< std::min(n,i+n/k); kk++)
 				B[j].push_back(kk);
 		}
 		else{
 			--j;
 			for(kk = i; kk< std::min(n,i+n/k); kk++)
 				B[j].push_back(kk);
 		}
 	}
}

bool write_results(const strstrMap &output, std::string file = "../output/results.csv"){
    std::stringstream ss;
    size_t j = 0;
    for(auto i=output.begin(); i != output.end();++i,++j )
    {
        std::cout<<"\""<<i->first<<"\" ,";
        if (i->first == "log") continue;
        if (j != output.size()-1)
            ss<< i->second<<",";  
        else
            ss<< i->second;
    }
    ss<<"\n";
    

    std::ofstream out(file.c_str(),std::ios::app);
    if(out.fail())
    {
        out.close();
        return false;
    }
    out << ss.str();
    out.close();
    return true;
}


size_t IO_construct_hypergraph(const char *hg_file, intintVec &hyperedges, intvec& init_nodes, size_t &maxid){
	/* Reads in the hypergraph file and initializes core values.
	 */
	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	fp = fopen(hg_file, "r");
	if (NULL == fp) {
		D(fprintf(stderr, "Cannot open hypergraph file!\n"));
		return 0;
	}

	char * ptr;
	size_t counter;
    maxid = 0;
	intset nodes;
	while ((read = getline(&line, &len, fp)) != -1) {
		ptr = strtok(line,",");
		counter = 0;

        // Read a hyperedge
        std::vector<size_t> hyperedge;
		while (ptr != NULL) {
			size_t nodeid = atol(ptr);
            hyperedge.push_back(nodeid);
			if (nodes.find(nodeid) == nodes.end()){
				init_nodes.push_back(nodeid);
				nodes.insert(nodeid);
			}
			maxid = std::max(nodeid,maxid);
			ptr = strtok(NULL,",");
			counter++;
		}
		// Add to e_id_toedge
		hyperedges.push_back(hyperedge);
	}
	if (line) {
		free(line);
	}
    return 1;
}
size_t init_cores(intintVec& hyperedges, intvec& min_e_hindex, intvec& llb, size_t& glb, intvec& init_nodes, intintMap& node_index, omp_lock_t lock[], size_t** nbrs_N,size_t** nbrs_F, size_t** inc_edges_N, size_t** inc_edges_F ) {
    /* Reads the incidence list init_nbr and initialises the data structure required for local core. */
    size_t N = init_nodes.size();
	double start_init = omp_get_wtime();
    for(size_t i = 0; i<N; i++) node_index[init_nodes[i]] = i; // initialize node_index
	double init1 = omp_get_wtime() - start_init;

    intvec nbrsizes(N); 
    std::unordered_map<size_t, intset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    std::vector<intvec> nbrs(N);
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
	size_t M = 0;
    std::vector<intvec> inc_edges(N); // i=node_id, value = vector of edge ids incident on node_id
	start_init = omp_get_wtime();
    for (size_t eid = 0; eid < hyperedges.size(); eid++){
		auto hype = hyperedges[eid];
        auto edge_sz = hype.size();
        sz_inc_edge += edge_sz;
		
		omp_init_lock(&(lock[eid]));
        size_t _min = INT_MAX;
		min_e_hindex[M] = _min; // initialize edge h_indices,
		M+=1;
        for(auto v_id: hype){
            auto j = node_index[v_id];  
            llb[j] = std::max(edge_sz - 1,llb[j]);
			// size_t i = node_index[v_id];
            inc_edges[j].push_back(eid);
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = intset();
                int _tmp_sz = 0;
                for (auto u: hype){
                    if (u!=v_id){
                        _tmp.insert(u);
                        _tmp_sz+=1;
                    }
                }
                init_nbr[v_id] = std::move(_tmp);
                nbrsizes[j] = _tmp_sz;
            }
            else{  // v_id exists in init_nbr map
                auto _tmp = &init_nbr[v_id];
                for (auto u: hype){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                nbrsizes[j] = init_nbr[v_id].size();
            }
        }
    }
	double init2 = omp_get_wtime() - start_init;

	start_init = omp_get_wtime();
	
	double init3 = omp_get_wtime() - start_init;

    for(auto i : init_nodes){
        sz_init_nbrs += init_nbr[i].size();
    }
    
	start_init = omp_get_wtime();
    *nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    (*nbrs_N)[0] = 0;
    *nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));

    glb = INT_MAX;
    int _i = 1, _index=0;
    for (auto node : init_nodes){
        (*nbrs_N)[_i] = (*nbrs_N)[_i-1] + init_nbr[node].size();
        _i++;
        size_t sz = 0;
        for (auto u: init_nbr[node]){
            // _tmp.push_back(node_index[u]);
            (*nbrs_F)[_index++] = node_index[u];
            sz+=1;
        }
        glb = std::min(glb, sz);
    }
	double initcsr1 = omp_get_wtime() - start_init;

	start_init = omp_get_wtime();
    // Calculate csr representation for incident edges
    *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    (*inc_edges_N)[0] = 0;
    *inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    _i = 1, _index=0;
    for(auto node : init_nodes){
        size_t j = node_index[node];
        (*inc_edges_N)[_i] = (*inc_edges_N)[_i-1] + inc_edges[j].size();
        _i++;
        for(size_t eid : inc_edges[j]){
            (*inc_edges_F)[_index++] = eid;
        }
    }
	double initcsr2 = omp_get_wtime() - start_init;
	std::cout<< "node_index: "<<init1<<" s\n";
	std::cout<< "nbrsize & nbr construction from Edges: "<<init2<<" s\n";
	std::cout << "incident edgeList construction: "<<init3<<" s\n";
	std::cout << "nbr_N & nbr_F: "<<initcsr1<<" s\n";
	std::cout<< "incedges_N & incedges_F: "<<initcsr2<<" s\n";
	return 1;
}


void write_core(std::vector<graph_node>& A, size_t N, const intvec& init_nodes, intintMap& node_index, intvec& node_index_index, std::string dataname = "output/hgcore"){
	std::string file = dataname+".core";
    std::stringstream ss;
    for(size_t i=0; i<N; i++)
    {
		if (A[i].valid) {
			auto u = init_nodes[i];
        	ss << u << "," << A[node_index_index[node_index[u]]].kcore << "\n";
		}
    }
    std::ofstream out(file.c_str());
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();
}

int main (int argc, char *argv[]) {
	if (argc < 2) {
		printf("usage: %s <graph-name> <num_threads>\n", argv[0]);
		return 1;
	}
	size_t working_threads = atoi(argv[2]);
	bool lbflag = true;
	if(argc>=3){
        std::string s = argv[3];
        if(s[0]=='0') lbflag = false;
    }

    strstrMap output;
    size_t maximum_id;
	intintVec hyperedges;
    auto hypergraph_file = dataset_to_filename[argv[1]];
	intvec init_nodes;
    if (IO_construct_hypergraph(hypergraph_file.c_str(), hyperedges, init_nodes,maximum_id) == 0){
        std::cout << "IO Error. Terminating..."<<"\n";
        return 1;
    }
	size_t N = init_nodes.size();
    //initialisation 
	double start_init = omp_get_wtime();
	std::vector<graph_node> A;	
	intvec min_e_hindex( hyperedges.size() );// min edge h-index for optimization II
	omp_lock_t *Elock = new omp_lock_t[hyperedges.size()]; // lock for hyperedges to apply optimization II without race-condition
	intvec llb(N,0);
	size_t glb;
	intintMap node_index;
	size_t* nbrs_N = NULL;
	size_t* nbrs_F = NULL;
	size_t* inc_edges_N = NULL;
	size_t* inc_edges_F = NULL;
	init_cores(hyperedges,min_e_hindex,llb,glb,init_nodes,node_index,Elock,&nbrs_N, &nbrs_F, &inc_edges_N, &inc_edges_F);
	// intintVec B;
	// intvec prefixsum_partition;
	// intvec node_index_index(node_index.size());
	// simple_loadbalance(A,N,working_threads,B,prefixsum_partition,node_index_index, hyperedges,node_index,nbrs_N);
	// if (lbflag)
	// 	re_order_loadbalance(N,working_threads,B,prefixsum_partition,node_index_index, hyperedges,node_index,nbrs_N);
	// else
	// 	no_loadbalance(N,working_threads,B,prefixsum_partition,node_index_index, hyperedges,node_index,nbrs_N);
	//  size_t index = 0;
	//  //std::cout<<"Allocating A:\n";
	//  for (size_t i = 0; i< working_threads; i++){
	//  	for(size_t j=0; j<B[i].size();j++){
	//  		graph_node new_node;
    // 		new_node.id = B[i][j];
	//  		new_node.kcore = (nbrs_N[B[i][j]+1]-nbrs_N[B[i][j]]);
	// 		//std::cout<<"x\n";
	//  		new_node.h_count.resize(new_node.kcore+1);
	// 		//std::cout<<"y\n";
	//  		new_node.num_nbrs = new_node.kcore;
	//  		new_node.active = 1;
	//  		new_node.valid = 1;
	// 		A.push_back(new_node);
	//  		node_index_index[B[i][j]] = index;
	//  		index++;
	//  	}
	//  	prefixsum_partition.push_back(index);
	// 	B[i].resize(0);	B[i].shrink_to_fit();
	// 	//std::cout<<"done A for Thread "<<i<<"\n";
	//  }
	//B.resize(0);
	//B.shrink_to_fit();
	double init_time = omp_get_wtime() - start_init;
    //core-computation
	// double core_start = omp_get_wtime();sta
	// size_t steps = compute_k_core(N, working_threads, A, node_index_index, min_e_hindex, Elock, prefixsum_partition, llb, glb, nbrs_N, nbrs_F, inc_edges_N,inc_edges_F,hyperedges, node_index, true);
	// double core_time = omp_get_wtime() - core_start;
    // printf("#Threads:%lu/Time:%f seconds/steps: %lu\n\n",working_threads, core_time,steps);
	printf("Init time: %lf\n",init_time);
	// write_core(A, N, init_nodes, node_index,  node_index_index, "../output/"+std::string(argv[1]));
	// if (lbflag)	output["algo"] = "LocalP(B+CSR)2";
	// else	output["algo"] = "LocalP(nolb)";
    // output["dataset"]=argv[1];
    // output["num_threads"] = std::to_string(working_threads);
    // output["execution time"]= std::to_string(core_time);
	// output["init_time"] = std::to_string(init_time);
    // output["total iteration"] = std::to_string(steps);
	// if (lbflag)	write_results(output,"../output/parout/results.csv");
	// else 	write_results(output,"../output/parout/results_nolb.csv");
	
	delete[] Elock;

	return 0;
}
