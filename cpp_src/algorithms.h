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
#include "hypergraph.h"

typedef  std::map<std::string, size_t> strIntMap;
typedef  std::map<std::string, std::vector<size_t>> strvIntMap;
typedef  std::map<std::string, std::set<size_t>> strsIntMap;
typedef  std::map<int, std::set<std::string>> intsStrMap;
typedef  std::map<std::string, std::vector<std::string>> strvStrMap;
typedef  std::map<size_t, std::vector<std::string>> intvStrMap;
typedef  std::vector<std::string> strvec;
typedef  std::set<std::string> strset;
typedef  std::vector<size_t> intvec;
typedef std::vector<std::pair<std::string, std::string>> strstrprvec;
typedef std::map<std::string, std::string> strstrMap;
typedef std::map <std::string,bool> strboolMap;


class Algorithm{
    Hypergraph hg;
    public:
    strIntMap core;
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
    void Peel(bool verbose = false);
    void EPeel(bool verbose = false);
    void local_core(bool log = false);
    void local_core_opt_core_correct(bool log = false);
    void local_core_omp(std::map<size_t, strvec >&, bool log = false);
    size_t iterative_core_correct(Hypergraph& H, std::string u, size_t core_u, strIntMap &hn);
    size_t iterative_core_correct_opt(Hypergraph& H, std::string u, size_t core_u);
    void printcore();
    bool write_results();
    void writecore();
    void writelog();
    void writeNbrQ();
};
void local_core( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log); 
void local_core_OPTI( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log);
void local_core_OPTII( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log);
void local_core_OPTIII( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log);
void local_core_OPTIV( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log);
void Peel( std::string dataset, std::map<size_t, strvec > e_id_to_edge, std::map<std::string, std::set<size_t> > inc_dict, strvec init_nodes, Algorithm& a);
void EPeel( std::string dataset, std::map<size_t, strvec > e_id_to_edge, std::map<std::string, std::set<size_t> > inc_dict, strvec init_nodes, Algorithm& a);



#endif