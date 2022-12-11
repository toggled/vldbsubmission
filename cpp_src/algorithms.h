#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include "hypergraph.h"
typedef  std::unordered_map<size_t, size_t> intIntMap;
typedef  std::map<std::string, size_t> strIntMap;
typedef  std::map<std::string, std::vector<size_t>> strvIntMap;
typedef  std::map<std::string, std::set<size_t>> strsIntMap;
typedef  std::map<size_t, std::set<std::string>> intsStrMap;
typedef  std::map<std::string, std::vector<std::string>> strvStrMap;
typedef  std::map<size_t, std::vector<std::string>> intvStrMap;
// typedef  std::unordered_map <size_t, std::vector<size_t>> uintvIntMap;
typedef  std::vector<std::string> strvec;
typedef  std::set<std::string> strset;
typedef  std::vector<size_t> intvec;
typedef std::vector<std::pair<std::string, std::string>> strstrprvec;
typedef std::map<std::string, std::string> strstrMap;
// typedef std::map <std::string,bool> strboolMap;
typedef std::unordered_map <size_t,bool> intboolMap;
typedef std::unordered_set<size_t> uintSet;
typedef std::vector<uintSet > uintsetvec;
typedef std::unordered_map<size_t, uintSet> intuSetintMap;
typedef std::map<std::string, std::string> strstrMap;
typedef std::vector< intvec > intintvec;
typedef std::pair<size_t,size_t> intpair;
class Algorithm{
    Hypergraph hg;
    public:
    intIntMap core;
    intIntMap secondcore;
    std::vector<std::vector<std::pair<size_t,size_t>>> score_m;
    double exec_time = 0;
    double core_exec_time = 0;
    double correction_time = 0;
    size_t nu_cu = 0;
    size_t num_nbr_queries = 0;
    strstrMap output;
    std::vector< strstrMap > hnlog;
    strstrMap timelogs;
    Algorithm( Hypergraph &H);
    ~Algorithm();
    // void Peel(bool verbose = false);
    // void EPeel(bool verbose = false);
    // void local_core(bool log = false);
    // void local_core_opt_core_correct(bool log = false);
    // void local_core_omp(std::map<size_t, strvec >&, bool log = false);
    // size_t iterative_core_correct(Hypergraph& H, std::string u, size_t core_u, strIntMap &hn);
    // size_t iterative_core_correct_opt(Hypergraph& H, std::string u, size_t core_u);
    void printcore();
    bool write_results();
    void writecore(std::string folder="../output/");
    void writelog();
    void writeNbrQ();
    void writekdcore(std::string folder="../output/");
};
// void local_core( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log); 
void print_bucket( intuSetintMap&, intvec&);
void local_core( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log); 
void local_core_OPTI( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void local_core_OPTII( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void local_core_OPTIII( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void local_core_OPTIV( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void Peel( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void EPeel( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void degreePeel( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void local_core_clique( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void local_core_bipartite( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void kdCorehybrid( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log);
#endif