#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include <algorithm>
#include <limits>
#include "hypergraph.h"

// typedef  std::map<std::string, size_t> strIntMap;
// typedef  std::map<std::string, std::vector<size_t>> strvIntMap;
// typedef  std::map<std::string, std::vector<std::string>> strvStrMap;
// typedef  std::map<size_t, std::vector<std::string>> intvStrMap;
// typedef  std::vector<std::string> strvec;
// typedef  std::vector<size_t> intvec;

// Hypergraph{

    // std::map<size_t, strvec > e_id_to_edge; // # key => hyperedge_id, value => List of vertices in a hyperedge
   

    // public:
    //  // Auxiliary variables
    // std::map<std::string, std::vector<size_t> > inc_dict;  //# key => node, value = List of incident hyperedge ids.
    // std::map<std::string, strvec > init_nbr;  //# key => node id, value => List of Neighbours.
    // strIntMap init_nbrsize; // # initial nbrhood sizes. 
    // strvec init_nodes;
    // /* data structures that make the algorithm more efficient */
    // std::map<size_t, size_t > edge_min_hindex;
    // strIntMap lub;
    // strIntMap llb;
    Hypergraph::Hypergraph(){}
    Hypergraph::~Hypergraph(){}

    void Hypergraph::addEdge(size_t id, strvec edge){
        e_id_to_edge[id] = std::move(edge);
    }
    void Hypergraph::initialise(){
        /*
        Initialise different variables e.g. number of neighbours, neighbour list etc.
        */
        for(auto elem: e_id_to_edge){
            auto edge_id = elem.first;
            // size_t edge_id = it->first;
            // for (auto u_it = it->second.cbegin(); u_it!= it->second.cend(); u_it++ ){
            //     auto v_id = *u_it;
            for(auto v_id: elem.second){
                // initalise init_nodes and incident dictionary
                if ( inc_dict.find(v_id) == inc_dict.end() ) { // v not find (first time)
                    inc_dict[v_id] = std::set<size_t>{edge_id};
                    init_nodes.push_back(v_id);
                }
                else 
                    inc_dict[v_id].insert(edge_id);
            }
        }
    }
    void Hypergraph::init_nbrs(){
        for(auto elem: e_id_to_edge){
            for(auto v_id: elem.second){
            // initialise number of neighbours and set of neighbours
                if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                    init_nbr[v_id] = strset();
                    for (auto u: elem.second){
                        if (u!=v_id){
                            init_nbr[v_id].insert(u);
                        }
                    }
                    init_nbrsize[v_id] =  init_nbr[v_id].size();
                }
                else{  // v_id exists in init_nbr map
                    for (auto u: elem.second){
                        if (u!=v_id){
                            init_nbr[v_id].insert(u);
                        }
                    }
                    init_nbrsize[v_id] = init_nbr[v_id].size();
                }
            }
        }   
    }
    void Hypergraph::init_edgehindex_lub(){
        /* Initialise data structures that make the algorithm more efficient */
        // std::map<size_t, size_t > edge_min_hindex; //# key = edge_id, value => min (h_index of vertices in hyperedge edge_id)
        for (auto v: init_nodes){
            size_t lub_v = lub[v];
            for (size_t e_id : inc_dict[v]){
                if (edge_min_hindex.find(e_id) == edge_min_hindex.end()){
                    // auto val = std::numeric_limits<size_t>::max();
                    edge_min_hindex[e_id] = lub_v;
                }
                else{
                    edge_min_hindex[e_id] = std::min(lub_v, edge_min_hindex[e_id]);
                }
            }
        }
    }
    void Hypergraph:: iterate_nbrs(std::string v, strvec & nbrs){
        /* Returns the set of neighbours of v.
            implements a traversal from vertex v to each of its neighbours in contrast to set in neighbors(). 
            It also returns an iterator. So it avoids creating the neighborhood list explicitely.
            Overall complexity: O(d(v) * |e_max|), where e_max = largest hyperedge 
       */
        auto incident_edges = inc_dict[v];  //# {O(1)}
        if (incident_edges.size()){
            strboolMap visited_dict;
            for (auto e_id : incident_edges){  //# { O(d(v)) }
                for (auto u : e_id_to_edge[e_id]){  //# { O(|e|)}
                    if (u != v){
                        if (visited_dict.find(u) == visited_dict.end()){ // u does not exists
                            visited_dict[u] = true;
                            nbrs.push_back(u);
                        }
                        if (!visited_dict[u]){
                            visited_dict[u] = true;
                            nbrs.push_back(u);
                        }
                    }
                }
            }
        }
    }
    void Hypergraph:: writeneighborhood(std::string file= "../output/log_dblpNv.csv"){
        std::stringstream ss;
        if(init_nbr.size()==0)  init_nbrs();
        for(auto pair: init_nbr){
            auto node = pair.first;
            ss<<node<<",";
            int _count = 0;
            int N = pair.second.size();
            for (auto nbr_v: pair.second){
                _count++;
                if(_count<N) ss<<nbr_v<<",";
                else    ss<<nbr_v<<"\n";
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
    size_t  Hypergraph:: get_number_of_nbrs(std::string v){
        std::set<std::string> nbrs;
        auto incident_edges = inc_dict[v];  //# {O(1)}
        if (incident_edges.size()){
            strboolMap visited_dict;
            for (auto e_id : incident_edges){  //# { O(d(v)) }
                for (auto u : e_id_to_edge[e_id]){  //# { O(|e|)}
                    if (u != v){
                        if (visited_dict.find(u) == visited_dict.end()){ // u does not exists
                            visited_dict[u] = true;
                            nbrs.insert(u);
                        }
                        if (!visited_dict[u]){
                            visited_dict[u] = true;
                            nbrs.insert(u);
                        }
                    }
                }
            }
        }
        return nbrs.size();
    }
    void Hypergraph::removeV_transform(std::string v, bool verbose){
        if (verbose)    std::cout << "remove v transform:\n";
        // For every edge e_id incident on v, remove e_id from every vertex u in edge[e_id] distinct from v
        for (auto e_id : inc_dict[v]){
            for(auto u: e_id_to_edge[e_id]){
                if (u==v)   continue;
                inc_dict[u].erase(e_id);
                if (verbose){
                    std::cout << "nbr -> "<<u<<"\n";
                    std::cout<<" rem: "<<e_id<<" from inc_dict["<<u<<"]\n";
                }
            }
        }
        inc_dict[v] = std::set<size_t>();
    }

    void Hypergraph::init_edgehindex_nbr(){
        /* Initialise edge minimum h index to be number of neighbours */
        for (auto v: init_nodes){
            size_t nbr_v = init_nbrsize[v];
            for (size_t e_id : inc_dict[v]){
                if (edge_min_hindex.find(e_id) == edge_min_hindex.end()){
                    // auto val = std::numeric_limits<size_t>::max();
                    edge_min_hindex[e_id] = nbr_v;
                }
                else{
                    edge_min_hindex[e_id] = std::min(nbr_v, edge_min_hindex[e_id]);
                }
            }
        }
    }

    void Hypergraph::compute_local_upperbound(){
        // # Computing global upper and lower bounds
        size_t glb = std::numeric_limits<size_t>::max();
        size_t gub = std::numeric_limits<size_t>::min();
        
        // # Auxiliary variables to assist computation of ub2
        strIntMap _inv_bucket;
        intvStrMap _bucket;
        
        // # Bucket initialisations
        for (auto v:init_nodes){
            auto len_neighbors_v = init_nbrsize[v];
            glb = std::min(glb, len_neighbors_v);
            gub = std::max(gub, len_neighbors_v);
            _inv_bucket[v] = len_neighbors_v;
            if (_bucket.find(len_neighbors_v) == _bucket.end())
                _bucket[len_neighbors_v] = strvec();
            _bucket[len_neighbors_v].push_back(v);
        }
        // # Week core-number computation
        for (auto k = glb; k < gub+1; k++){
            while (_bucket[k].size() != 0){
                auto v = _bucket[k].back();
                _bucket[k].pop_back();
                lub[v] = k;
                for (auto u : init_nbr[v]){
                    if (lub.find(u) == lub.end()){
                        auto max_value = std::max(_inv_bucket[u] - 1, k);
                        auto u_it = std::find(_bucket[_inv_bucket[u]].begin(),_bucket[_inv_bucket[u]].end(),u);
                        _bucket[_inv_bucket[u]].erase(u_it);
                        if(_bucket.find(max_value) == _bucket.end()){
                             _bucket[max_value] = strvec();
                        }    
                        _bucket[max_value].push_back(u);
                        _inv_bucket[u] = max_value;
                    }
                }
            }
        }

    }
    
    void Hypergraph::compute_local_lowerbound(){
        /*
        Computes local lower bound for every node
        */

        glb = std::numeric_limits<size_t>::max();
        gub = std::numeric_limits<size_t>::min();
        // # Computing global lower bounds
        for (auto v :init_nodes){
            auto len_neighbors_v = init_nbrsize[v];
            glb = std::min(glb,len_neighbors_v);
            gub = std::max(gub, len_neighbors_v);
        }
        // # Computing Local lower bound 
        for (auto v :init_nodes){
            auto _maxv = std::numeric_limits<size_t>::min();
            for (auto e_id : inc_dict[v]){
                size_t vec_sz = e_id_to_edge[e_id].size();
                _maxv = std::max(_maxv, vec_sz - 1);
            }
            llb[v] = std::max(_maxv, glb);
        }
    }

    size_t Hypergraph::get_min_hindex(size_t e_id){
        return edge_min_hindex[e_id];
    }

    void Hypergraph::update_min_hindex(std::string v, size_t h_v){
        /* Given a new value of h_v, update h_min for its incident hyperedges, whenever appropriate */
        for (auto e_id : inc_dict[v]){
            if (get_min_hindex(e_id) > h_v)
                edge_min_hindex[e_id] = h_v;
        }
    }

    void Hypergraph::test_initialisation(){
        // print init_nodes
        for (auto v: init_nodes) std::cout<< v<<" ";
        std::cout<<"\n";
        // print incident edges on v for all v
        for(auto it = inc_dict.cbegin(); it != inc_dict.cend(); ++it){
            std::cout<< it->first<<": ";
            for (auto e_id: it->second){
                std::cout<<e_id<<" ";
            }
            std::cout<<'\n';
        }
        // print nbr(v) => len(nbr(v)) for all v
        for(auto it = init_nbr.cbegin(); it != init_nbr.cend(); ++it){
            std::cout<< it->first<<": ";
            for (auto nbr_v: it->second){
                std::cout<<nbr_v<<" ";
            }
            std::cout << "=>"<<init_nbrsize[it->first];
            std::cout<<'\n';
        }

        std::cout << "Edge min hindex: \n";
        for(const auto& elem : edge_min_hindex)
        {
        std::cout << elem.first << "->"<<elem.second<<"\n";
        }

        std::cout << "local upper bounds: \n";
        for(const auto& elem : lub)
        {
        std::cout << elem.first << "->"<<elem.second<<"\n";
        }

        std::cout << "local lower bounds: \n";
        for(const auto& elem : llb)
        {
        std::cout << elem.first << "->"<<elem.second<<"\n";
        }
    }

    void Hypergraph::printHypergraph(){
        /* 
        Prints the edge list 
        */
        // for(auto it = e_id_to_edge.cbegin(); it != e_id_to_edge.cend(); ++it){
        //     std::cout << it->first << ":";
        //     for (auto u_it = it->second.cbegin(); u_it!= it->second.cend(); u_it++ ){
        //         std::cout << *u_it<<" ";
        //     }   
        //     std::cout<<"\n";
        // }
        std::cout<< "hyperedges: \n";
        for (auto elem: e_id_to_edge){
            std::cout <<elem.first<<":";
            for(auto u: elem.second){
                std::cout << u <<" ";
            }
            std::cout<<"\n";
        }
        std::cout<<" incidence: \n";
        for(auto it = inc_dict.cbegin(); it != inc_dict.cend(); ++it){
            std::cout<< it->first<<": ";
            for (auto e_id: it->second){
                std::cout<<e_id<<" ";
            }
            std::cout<<'\n';
        }
    }
// };