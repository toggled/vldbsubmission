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
#include <unordered_set>
#include <array>
typedef  std::vector<std::string> strvec;
typedef  std::vector<size_t> intvec;
typedef  std::set<size_t> intset;
typedef	 std::vector<intset> vintset;
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
            {"pref", "../data/datasets/synthetic/pref_1000000_3_1.hyp"},
			{"default" , "../data/datasets/synthetic/default.hyp"}
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

bool write_results(const strstrMap &output, std::string file = "../output/inresults.csv"){
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
    std::cout<<"\n";

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

bool LCCSAT_csr(const size_t inc_edges_F[], const size_t inc_edges_N[],intvec& min_hindices, const intintVec& edges, size_t u_id, size_t core_u){
	std::set <size_t> Nplus;
    for (size_t i = inc_edges_N[u_id]; i<inc_edges_N[u_id+1]; i++){
        size_t e_id = inc_edges_F[i];
        if (min_hindices[e_id] >= core_u){
            for(auto v_id: edges[e_id]){
                Nplus.insert(v_id);
            }
            if(Nplus.size()-1>=core_u)	return true;
        }
    }
    auto sz_Nplus = Nplus.size();
    if (sz_Nplus>0){
        if (sz_Nplus-1 >= core_u)
            return true;
        else
            return false;
    }
    else
        return false;
}

bool LCCSAT_opt(std::vector<intvec> & min_hindex_to_edge, const intintVec& edges, size_t core_u, std::set<size_t> & Nplus){
    for(auto eid : min_hindex_to_edge[core_u]){
        for(auto v : edges[eid]){
            Nplus.insert(v);
        }
    }
    auto sz_Nplus = Nplus.size();
    if (sz_Nplus>0){
        if (sz_Nplus-1 >= core_u)
            return true;
        else
            return false;
    }
    else
        return false;
}

size_t core_correct_csr(const size_t inc_edges_F[], const size_t inc_edges_N[],intvec& min_hindices, const std::vector<intvec>&edges, size_t u_id, size_t core_u){  
        core_u = core_u - 1;
        std::vector<intvec> min_hindex_to_edge(core_u+1);
        for (size_t i = inc_edges_N[u_id]; i<inc_edges_N[u_id+1]; i++){
            size_t eid = inc_edges_F[i];
            if(min_hindices[eid]>=core_u){
                min_hindex_to_edge[core_u].push_back(eid);
            }else 
                {
                min_hindex_to_edge[min_hindices[eid]].push_back(eid);
            }
        }

        std::set<size_t> Nplus;
        
        while (LCCSAT_opt(min_hindex_to_edge, edges, core_u, Nplus) == false) {
            core_u = core_u - 1;
        }

        return core_u;
}

size_t process_message(size_t id, std::vector<graph_node>& A) {
	/* computes h-index.
	 */

	size_t ret = 0;
	A[id].active = 0;

	// ----- csr + remapped to [0,n<<N]  ------
	size_t t = 0;
	if (A[id].num_nbrs!=0){
		size_t paper = 0;
		size_t n = A[id].num_nbrs;
        for(size_t i = n; i >= 0; --i){
            paper += A[id].h_count[i];
            if(paper >= i){
                t = i;
				break;
			}
        }
	}
	if (t<A[id].kcore){
		ret = 1;
		A[id].active = 1;
		D(printf("%d is setting core from %d to %d\n", id, A[id].kcore, t));
		A[id].kcore = t;
	}
	else{
		D(printf("%d not active\n",id));
	}
	return ret;
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
size_t init_cores(intintVec& hyperedges, intvec& min_e_hindex, intvec& llb, size_t& glb, intvec& init_nodes, intintMap& node_index, omp_lock_t lock[], size_t** nbrs_N,size_t** nbrs_F, size_t** inc_edges_N, size_t** inc_edges_F, size_t working_threads, omp_lock_t Vlock[]) {
    /* Reads the incidence list init_nbr and initialises the data structure required for local core. */
    bool log = false;
	size_t N = init_nodes.size();
	// double start_init = omp_get_wtime();
	// double init1 = omp_get_wtime() - start_init;
	// for(size_t i = 0; i<N; i++) {
	// 	std::cout<< i<<": init_nodes[i] = "<<init_nodes[i]<<" mapped to node_index: "<<node_index[init_nodes[i]]<<"\n";
	// }
    // intvec nbrsizes(N); 
    // std::unordered_map<size_t, intset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
	double start_init = omp_get_wtime();
	// vintset init_nbr(N, intset({}));
	std::vector< std::unordered_set<size_t>> init_nbr(N, std::unordered_set<size_t>{});
	// std::vector< std::unordered_map<size_t,size_t>> init_nbr(N);
	// std::vector<intvec> nbrs(N);
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
	size_t M = 0;
    std::vector<intvec> inc_edges(N); // i=node_id, value = vector of edge ids incident on node_id
	// std::vector<bool> traversed(N,false);
    // for (size_t eid = 0; eid < hyperedges.size(); eid++){
	// 	auto hype = hyperedges[eid];
    //     auto edge_sz = hype.size();
    //     sz_inc_edge += edge_sz;
		
	// 	omp_init_lock(&(lock[eid]));
    //     size_t _min = INT_MAX;
	// 	min_e_hindex[M] = _min; // initialize edge h_indices,
	// 	M+=1;
    //     for(auto v_id: hype){
    //         auto j = node_index[v_id];  
    //         llb[j] = std::max(edge_sz - 1,llb[j]);
	// 		// size_t i = node_index[v_id];
    //         inc_edges[j].push_back(eid);
    //         // if (!traversed[j]) { // first insertion of v_id to init_nbr map
    //         //     auto _tmp = intset();
    //         //     int _tmp_sz = 0;
    //         //     for (auto u: hype){
    //         //         if (u!=v_id){
    //         //             _tmp.insert(u);
    //         //             _tmp_sz+=1;
    //         //         }
    //         //     }
    //         //     // init_nbr[v_id] = std::move(_tmp);
	// 		// 	init_nbr[j] = std::move(_tmp);
    //         //     nbrsizes[j] = _tmp_sz;
	// 		// 	traversed[j] = true;
    //         // }
    //         // else{  // v_id exists in init_nbr map
    //             // auto _tmp = &init_nbr[v_id];
	// 			auto _tmp = &init_nbr[j];
    //             for (auto u: hype){
    //                 if (u!=v_id){
    //                     _tmp->insert(u);
    //                 }
    //             }
    //             // nbrsizes[j] = init_nbr[v_id].size();
	// 			// nbrsizes[j] = init_nbr[j].size();
    //         // }
    //     }
    // }
	size_t eid;
	#pragma omp parallel for schedule(dynamic) num_threads(working_threads)
	for (eid = 0; eid < hyperedges.size(); eid++){
		auto hype = hyperedges[eid];
                auto edge_sz = hype.size();
                // sz_inc_edge += edge_sz;
		min_e_hindex[eid] = INT_MAX;
		for(auto v_id: hype){
            auto j = node_index[v_id]; 
			omp_set_lock(&(Vlock[j]));     
            		llb[j] = std::max(edge_sz - 1,llb[j]);
			inc_edges[j].push_back(eid);
			auto _tmp = &init_nbr[j];
			for (auto u: hype){
				if (u!=v_id){
					_tmp->insert(u);
				}
				// if (u!=v_id){
				// 	// not found
				// 	if (_tmp->find(u) == _tmp->end()){
				// 		(*_tmp)[u] = 1;
				// 	}
				// }
			}
			omp_unset_lock(&(Vlock[j])); 
		}
	}
	double init2 = omp_get_wtime() - start_init;
	// std::cout<<"init_nbr: \n";
	// for(int j=0; j<N; j++){
	// 	std::cout<<j<<": ";
	// 	for(auto u: init_nbr[j])	std::cout<<u<<" ";
	// 	std::cout<<"\n";
	// }	
	// std::cout<<"inc_edges: \n";
	// for(int j=0; j<N; j++){
	// 	std::cout<<j<<": ";
	// 	for(auto eid: inc_edges[j])	std::cout<<eid<<" ";
	// 	std::cout<<"\n";
	// }	
	// start_init = omp_get_wtime();
	
	// double init3 = omp_get_wtime() - start_init;

	glb = INT_MAX;
    // for(auto i : init_nodes){
	// 	auto sz = init_nbr[i].size();
    //     sz_init_nbrs += sz;
	// 	glb = std::min(glb, sz);
    // }
	start_init = omp_get_wtime();
	for(int _i = 0; _i< N; _i ++){
		auto sz = init_nbr[_i].size();
        sz_init_nbrs += sz;
		glb = std::min(glb, sz);
   	 }
	// for (auto hyp: hyperedges){
	for (eid = 0; eid < hyperedges.size(); eid++){
		auto hyp = hyperedges[eid];
		sz_inc_edge += hyp.size();
		omp_init_lock(&(lock[eid]));
	}
	double prefixsum_tm = omp_get_wtime()-start_init;
    
	start_init = omp_get_wtime();
    *nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    (*nbrs_N)[0] = 0;
    *nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
	*inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    (*inc_edges_N)[0] = 0;
    *inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
//    _index=0;
	int _i;
	if (log){
		std::cout<<"init_nodes: \n";
		for(int i = 0; i<N; i++) 	std::cout<<init_nodes[i]<<" ";
		std::cout<<"\n----\n";
		std::cout<<"node_index: \n";
		for(int i = 0; i<N; i++) 	std::cout<<node_index[init_nodes[i]]<<" ";
		std::cout<<"\n----\n";
	}
    for (_i = 1; _i<= N; _i ++){
		// auto node = init_nodes[_i-1];
		// std::cout<<node<<" "<<_i<<"/"<<N<<": "<<omp_get_thread_num()<<"\n";
        (*nbrs_N)[_i] = (*nbrs_N)[_i-1] + init_nbr[_i-1].size();
		(*inc_edges_N)[_i] = (*inc_edges_N)[_i-1] + inc_edges[_i-1].size();
		// std::cout<<node<<" "<<(*nbrs_N)[_i] <<" "<<_i<<"/"<<N<<": "<<omp_get_thread_num()<<"\n";
    }
	#pragma omp parallel for schedule(dynamic) num_threads(working_threads)
	for (int _i = 1; _i< N; _i ++){
		// auto node = init_nodes[_i-1];
		// auto it = init_nbr[node].begin();
				// std::cout<<node<<" "<<(*nbrs_N)[_i] <<" it "<<&it<<" "<<_i<<"/"<<N<<": "<<omp_get_thread_num()<<"\n";

		// _tmp.push_back(node_index[u]);
		// std::cout<<(*nbrs_N)[_i-1]<<" - "<<(*nbrs_N)[_i];
		// for(int _index = (*nbrs_N)[_i-1]; _index < (*nbrs_N)[_i]; _index++){
		// 	// std::cout<<"["<<*it<<" , "<<node_index[*it]<<"]";
		// 	(*nbrs_F)[_index] = node_index[*it];
		// 	it++;
		// }
		int _index = (*nbrs_N)[_i-1];
		// for(auto u: init_nbr[node]){
		for(auto u: init_nbr[_i-1]){
			(*nbrs_F)[_index++] = node_index[u];
			// (*nbrs_F)[_index++] = node_index[u.first];
		}
		_index = (*inc_edges_N)[_i-1];
		for(size_t eid : inc_edges[_i-1])
			(*inc_edges_F)[_index++] = eid;
		// std::cout<<"\n";
	}
	// log = true;
	double initcsr1 = omp_get_wtime() - start_init;
	if (log){
		std::cout<<"nbrs_F\n";
		for(int i = 0; i<20; i++)	std::cout<<(*nbrs_F)[i]<<" ";
		std::cout<<"\n";
		std::cout<<"init_nbr: \n";
		for (int i = 0; i< N; i++){
			auto node = init_nodes[i];
			std::cout<< node<<": ";
			// for(auto u: init_nbr[i])	std::cout<<u<<" ";
			// std::cout<<"\n";
		}
		std::cout<<"---\n";
	}
	start_init = omp_get_wtime();
    // Calculate csr representation for incident edges
    // *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    // (*inc_edges_N)[0] = 0;
    // *inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    // _i = 1;
	// int _index=0;
    // for(auto node : init_nodes){
    //     size_t j = node_index[node];
    //     (*inc_edges_N)[_i] = (*inc_edges_N)[_i-1] + inc_edges[j].size();
    //     _i++;
    //     for(size_t eid : inc_edges[j]){
    //         (*inc_edges_F)[_index++] = eid;
    //     }
    // }
	// int _j;
	// for (_j = 1; _j<= N; _j++)
	// 	(*inc_edges_N)[_j] = (*inc_edges_N)[_j-1] + inc_edges[_j-1].size();
	
	// #pragma omp parallel for schedule(dynamic) num_threads(working_threads)
	// for (_j = 1; _j<= N; _j++){
	// 	// (*inc_edges_N)[_j] = (*inc_edges_N)[_j-1] + inc_edges[_j-1].size();
	// 	int _index = (*inc_edges_N)[_j-1];
	// 	for(size_t eid : inc_edges[_j-1])
	// 		(*inc_edges_F)[_index++] = eid;
	// }
	double initcsr2 = omp_get_wtime() - start_init;
	
	// std::cout<<"incedges_N & incedges_F: \n";
	// for(int i = 0; i<N+1; i++) std::cout<< (*inc_edges_N)[i]<<" ";
	// std::cout<<"\n";
	// for(int i = 0; i<sz_inc_edge; i++) std::cout<< (*inc_edges_F)[i]<<" ";
	// std::cout<<"\n";
	// std::cout<<"nbrs_N & nbrs_F: \n";
	// for(int i = 0; i<N+1; i++) std::cout<< (*nbrs_N)[i]<<" ";
	// std::cout<<"\n";
	// for(int i = 0; i<sz_init_nbrs; i++) std::cout<< (*nbrs_F)[i]<<" ";
	// std::cout<<"\n";

	// std::cout<< "node_index: "<<init1<<" s\n";
	std::cout<< "nbrsize & nbr construction from Edges: "<<init2<<" s\n";
	//std::cout << "incident edgeList construction: "<<init3<<" s\n";
	std::cout<< "prefix sum time: "<<prefixsum_tm<<" s\n";
	std::cout << "nbr_N & nbr_F: "<<initcsr1<<" s\n";
	std::cout<< "incedges_N & incedges_F: "<<initcsr2<<" s\n";
	return 1;
}
int compute_k_core(size_t n, size_t working_threads, std::vector<graph_node>& A, const intvec& node_index_index, intvec& min_e_hindex, omp_lock_t elock[], const intvec& prefixsum_partition, const intvec& llb, int glb, const size_t nbrs_N[], const size_t nbrs_F[], const size_t inc_edges_N[], const size_t inc_edges_F[], const intintVec& hyperedges, intintMap &node_index, bool core_cct = false){
	// std::cout<<"compute_k_core\n";
	size_t supersteps = 0;
	if (working_threads > n) {
		// The chances of this happening are very slim.
		working_threads = 1;
	}
	D(printf("Initializing with %d\n", working_threads));

	size_t offsets[working_threads];
	size_t i;
	for (i = 0; i < working_threads; i++) {
		offsets[i] = 0;
	}
	
	size_t continue_itr = 1;
	while (continue_itr) {
		D(std::cout <<"Superstep: "<<supersteps<<"\n");
		#pragma omp parallel num_threads(working_threads)
		{
			size_t tid = omp_get_thread_num();
			size_t slice_start,slice_end;
			if (tid==0)	{
				slice_start = 0;
				slice_end = prefixsum_partition[tid];
			}
			else{ 
				slice_start = prefixsum_partition[tid-1];
				slice_end = prefixsum_partition[tid];
			}
			D(std::cout<<"Thread: "<<tid<<": "<<slice_start<<" "<<slice_end<<"\n");
			// #pragma omp for schedule(dynamic) // Comment  for not-dyn sched.
			for (size_t z = slice_start; z<slice_end; z++) {
				/* Want to stride across assigned blocks.
				 */
				if (A[z].valid) {
					if (0 == supersteps) {
						std::fill(A[z].h_count.begin(), A[z].h_count.end(), 0);
						if(core_cct){
							// updates min-hindex[e] to min(current min-hindex[e],A[z].kcore) for every e incident on z
							size_t l = inc_edges_N[A[z].id];
							size_t h = inc_edges_N[A[z].id+1];
							for(size_t i = l; i< h; i++){
								size_t eid = inc_edges_F[i];
								omp_set_lock(&(elock[eid]));
								min_e_hindex[eid] = std::min(A[z].kcore,min_e_hindex[eid]);
								omp_unset_lock(&(elock[eid]));
							}
						}
					} 
					else 
					{
						std::fill(A[z].h_count.begin(), A[z].h_count.end(), 0);
						size_t l = nbrs_N[A[z].id];
						size_t h = nbrs_N[A[z].id+1];
						for (size_t curr = l;  curr< h; curr++) {
								size_t nbr = node_index_index[nbrs_F[curr]];
								size_t inc_pos = std::min(A[z].num_nbrs, A[nbr].kcore);
								A[z].h_count[inc_pos]++;
						}
						offsets[tid] += process_message(z, A);
						if(core_cct){
							if (A[z].active){
								size_t l = inc_edges_N[A[z].id];
								size_t h = inc_edges_N[A[z].id+1];
								for(size_t i = l; i< h; i++){
									size_t eid = inc_edges_F[i];
									omp_set_lock(&(elock[eid]));
									min_e_hindex[eid] = std::min(A[z].kcore,min_e_hindex[eid]);
									omp_unset_lock(&(elock[eid]));
								}
							}
							bool lccsat = LCCSAT_csr(inc_edges_F,inc_edges_N,min_e_hindex,hyperedges, A[z].id, A[z].kcore);
							// std::cout <<"LCCSAT("<<z<<"/"<<A[z].id<<"):"<<lccsat<<"\n";
							if (lccsat == false){
								offsets[tid] += 1;
								auto hhatn = core_correct_csr(inc_edges_F,inc_edges_N,min_e_hindex,hyperedges,A[z].id,A[z].kcore);
								D(printf("%d: core corrected from %d to %d\n",A[z].id, A[z].kcore,hhatn));
								A[z].kcore  = hhatn;
								A[z].active = 1;
								size_t l = inc_edges_N[A[z].id];
								size_t h = inc_edges_N[A[z].id+1];
								for(size_t i = l; i< h; i++){
									size_t eid = inc_edges_F[i];
									omp_set_lock(&(elock[eid]));
									min_e_hindex[eid] = std::min(A[z].kcore,min_e_hindex[eid]);
									omp_unset_lock(&(elock[eid]));
								}
							}
						}
					}
				}
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			D(printf("continue/not-continue decision\n"));
			if (supersteps != 0) {
				continue_itr = 0;
				for (i = 0; i < working_threads; i++) {
					if (offsets[i] != 0) {
						D(printf("Need to continue\n"));
						continue_itr = 1;
					}
					offsets[i] = 0;
				}
			}

			supersteps++;
		}
	}
	return supersteps;
}
void write_core(std::vector<graph_node>& A, size_t N, const intvec& init_nodes, intintMap& node_index, intvec& node_index_index, std::string dataname = "output/hgcore"){
	std::string file = dataname+"_inopt.core";
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
	// std::cout<<hypergraph_file<<"\n";
	intvec init_nodes;
	intintMap node_index;
    if (IO_construct_hypergraph(hypergraph_file.c_str(), hyperedges, init_nodes,maximum_id) == 0){
        std::cout << "IO Error. Terminating..."<<"\n";
        return 1;
    }
	size_t N = init_nodes.size();
	std::vector<graph_node> A;	
	intvec min_e_hindex( hyperedges.size() );// min edge h-index for optimization II
	omp_lock_t *Elock = new omp_lock_t[hyperedges.size()]; // lock for hyperedges to apply optimization II without race-condition
	omp_lock_t *Vlock = new omp_lock_t[N];
	intvec llb(N,0);
	size_t glb;
	for(size_t i = 0; i<N; i++) {
		node_index[init_nodes[i]] = i; // initialize node_index
		omp_init_lock(&(Vlock[i]));
	}
	// for(int i = 0; i<N; i++) std::cout<< init_nodes[i]<<"\n";
    //initialisation 
	double start_init = omp_get_wtime();
	size_t* nbrs_N = NULL;
	size_t* nbrs_F = NULL;
	size_t* inc_edges_N = NULL;
	size_t* inc_edges_F = NULL;
	double start_initcores = omp_get_wtime();
	init_cores(hyperedges,min_e_hindex,llb,glb,init_nodes,node_index,Elock,&nbrs_N, &nbrs_F, &inc_edges_N, &inc_edges_F, working_threads, Vlock);
	double initctime = omp_get_wtime() - start_initcores;
	intintVec B;
	intvec prefixsum_partition;
	intvec node_index_index(node_index.size());
	double start_lb = omp_get_wtime();
	// simple_loadbalance(A,N,working_threads,B,prefixsum_partition,node_index_index, hyperedges,node_index,nbrs_N);
	if (lbflag)
		re_order_loadbalance(N,working_threads,B,prefixsum_partition,node_index_index, hyperedges,node_index,nbrs_N);
	else
		no_loadbalance(N,working_threads,B,prefixsum_partition,node_index_index, hyperedges,node_index,nbrs_N);
	double lbtime = omp_get_wtime() - start_lb;
	double start_all = omp_get_wtime();
	size_t index = 0;
	//  std::cout<<"Allocating A:\n";
	for (size_t i = 0; i< working_threads; i++){
	 	for(size_t j=0; j<B[i].size();j++){
	 		graph_node new_node;
    		new_node.id = B[i][j];
	 		new_node.kcore = (nbrs_N[B[i][j]+1]-nbrs_N[B[i][j]]);
			//std::cout<<"x\n";
	 		new_node.h_count.resize(new_node.kcore+1);
			//std::cout<<"y\n";
	 		new_node.num_nbrs = new_node.kcore;
	 		new_node.active = 1;
	 		new_node.valid = 1;
			A.push_back(new_node);
	 		node_index_index[B[i][j]] = index;
	 		index++;
	 	}
	 	prefixsum_partition.push_back(index);
		B[i].resize(0);	B[i].shrink_to_fit();
		// std::cout<<"done A for Thread "<<i<<"\n";
	 }
	// for(auto p: prefixsum_partition)	std::cout<<p<<"\n";
	B.resize(0);
	B.shrink_to_fit();
	double arrayofstructtime = omp_get_wtime() - start_all;
	double init_time = omp_get_wtime() - start_init;
    //core-computation
	double core_start = omp_get_wtime();
	//size_t steps = compute_k_core(N, working_threads, A, node_index_index, min_e_hindex, Elock, prefixsum_partition, llb, glb, nbrs_N, nbrs_F, inc_edges_N,inc_edges_F,hyperedges, node_index, true);
	double core_time = omp_get_wtime() - core_start;
    // printf("#Threads:%lu/Time:%f seconds/steps: %lu\n\n",working_threads, core_time,steps);
	printf("Init time: %lf\n",init_time);
	printf("initcore(): %lf\n",initctime);
	printf("arrayofst time: %lf\n",arrayofstructtime);
	printf("lb time: %lf\n",lbtime);
	write_core(A, N, init_nodes, node_index, node_index_index, "../output/parout/"+std::string(argv[1])+"_"+argv[2]);
	if (lbflag)	output["algo"] = "LocalP(B+CSR)2";
	else	output["algo"] = "LocalP(nolb)";
    output["dataset"]=argv[1];
    output["num_threads"] = std::to_string(working_threads);
    output["execution time"]= std::to_string(core_time);
	output["init_time"] = std::to_string(init_time);
    // output["total iteration"] = std::to_string(steps);
	if (lbflag)	write_results(output,"../output/parout/inresults.csv");
	else 	write_results(output,"../output/parout/inresults_nolb.csv");
	
	delete[] Elock;
	delete[] Vlock;
	return 0;
}
