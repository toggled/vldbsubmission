#ifndef HYPERGRAPH_H
#define HYPERGRAPH_H
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <string> 
#include <set>
#include <algorithm>
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>

typedef std::unordered_map<std::string, size_t> strInthashMap;
typedef  std::map<std::string, size_t> strIntMap;
typedef std::map<size_t, size_t > intIntMap;
typedef  std::map<std::string, std::vector<size_t>> strvIntMap;
typedef  std::map<std::string, std::vector<std::string>> strvStrMap;
typedef  std::map<size_t, std::vector<std::string>> intvStrMap;
typedef  std::set<std::string> strset;
typedef  std::vector<std::string> strvec;
typedef  std::vector<size_t> intvec;
typedef std::map<std::string, std::string> strstrMap;
typedef std::map <std::string,bool> strboolMap;
class Hypergraph{
    public:
    std::string dataset;
     // Auxiliary variables
    std::map<size_t, strvec > e_id_to_edge; // # key => hyperedge_id, value => List of vertices in a hyperedge
    std::map<std::string, std::set<size_t> > inc_dict;  //# key => node, value = List of incident hyperedge ids.
    std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    strIntMap init_nbrsize; // # initial nbrhood sizes. 
    strvec init_nodes;
    /* data structures that make the algorithm more efficient */
    std::map<size_t, size_t > edge_min_hindex;
    strIntMap lub;
    strIntMap llb;
    size_t glb;
    size_t gub;
    Hypergraph();
    ~Hypergraph();

    void addEdge(size_t id, strvec edge);
    void initialise();
    void init_nbrs();
    void init_edgehindex_lub();
    void init_edgehindex_nbr();
    void iterate_nbrs(std::string, strvec & );
    size_t get_number_of_nbrs(std::string );
    void removeV_transform(std::string, bool);

    void compute_local_upperbound();
    void compute_local_lowerbound();
    size_t get_min_hindex(size_t e_id);
    void update_min_hindex(std::string v, size_t h_v);


    void test_initialisation();
    void printHypergraph();
    void writeneighborhood(std::string);
};
#endif