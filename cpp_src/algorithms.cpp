#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include <algorithm>
#include "hypergraph.h"
#include "algorithms.h"
#include "utils.h"

//--------------------------------------------------------------------- Utility functions ------------------------------------------------------------------------------

void Algorithm::printcore(){
    std::cout << "core: \n";
    for(const auto& elem : core)
    {
    std::cout << elem.first << "->"<<elem.second<<"\n";
    }
}
void Algorithm::writecore(){
    // std::cout << "core: \n";
    std::string file = "../output/core_"+output["algo"]+"_"+hg.dataset+".csv";
    std::stringstream ss;
    for(auto elem: core)
    {
        ss << elem.first << "," << elem.second << "\n";
    }
    std::ofstream out(file.c_str());
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();
}
bool Algorithm::write_results(){
    std::string file = "../output/results.csv";
    std::stringstream ss;
    
    // // print the order in which the output keys are saved in the csv file.
    // for(auto i=output.begin(); i != output.end();++i )
    //     std::cout<< i->first<<",";
    
    size_t j = 0;
    for(auto i=output.begin(); i != output.end();++i,++j )
    {
        // std::cout<<"\""<<i->first<<"\" ,";
        if (i->first == "log") continue;
        // ss << i->first << "," << i->second << "\n";
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

void Algorithm::writeNbrQ(){
    std::string file = "../output/nbrq_results.csv";
    std::stringstream ss;
    ss<<output["algo"]<<" "<<hg.dataset<<" "<<output["execution time"]<<" "<<num_nbr_queries<<"\n";
    std::ofstream out(file.c_str(),std::ios::app);
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();
}

void  Algorithm::writelog(){
    std::string file = "../output/log_"+output["algo"]+"_"+hg.dataset+".csv";
    std::stringstream ss;
    
    size_t i = 0;
    for (auto hn_i:hnlog){
        size_t j = 0;
        // writing column names (only once)
        if (i==0){
            for(auto key:hn_i){
                if (j != hn_i.size()-1)
                    ss<<key.first<<",";
                else
                    ss<<key.first;
                j++;
            }
            ss<<"\n";
        }
        //writing values (at each iteration)
        j = 0;
        for(auto key:hn_i){
            if (j != hn_i.size()-1)
                ss<<key.second<<",";
            else
                ss<<key.second;
            j++;
         }
         ss<<"\n";
         i++;
    }

    std::ofstream out(file.c_str());
    if(out.fail())
    {
        out.close();
        std::cout<<"log writing failed\n";
        // return false;
    }
    out << ss.str();
    out.close();
    // return true;
}
Algorithm::Algorithm(Hypergraph &H){
    hg = H;
    output["dataset"] = H.dataset;
}
Algorithm::~Algorithm(){}

// ------------------------------------------------------------------------ Peel ------------------------------------------------------------------------------------

void iterate_nbrs(std::string v, strvec & nbrs, std::map<std::string, std::set<size_t> > & inc_dict, std::map<size_t, strvec > &e_id_to_edge){
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

void removeV_transform(std::string v, std::map<std::string, std::set<size_t> > & inc_dict, std::map<size_t, strvec > &e_id_to_edge){
    // For every edge e_id incident on v, remove e_id from every vertex u in edge[e_id] distinct from v
    for (auto e_id : inc_dict[v]){
        for(auto u: e_id_to_edge[e_id]){
            if (u==v)   continue;
            inc_dict[u].erase(e_id);
        }
    }
    inc_dict[v] = std::set<size_t>();
}

size_t get_number_of_nbrs(std::string v, std::map<std::string, std::set<size_t> > & inc_dict, std::map<size_t, strvec > &e_id_to_edge){
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

void Peel(std::string dataset, std::map<size_t, strvec > e_id_to_edge, std::map<std::string, std::set<size_t> > inc_dict, strvec init_nodes, Algorithm& a, bool log){
    a.output["algo"] = "Peel";
    clock_t start, end;
    // strIntMap nbrquery_stat;
    start = clock();
    intsStrMap bucket;
    strIntMap inverse_bucket;

    clock_t s_tm; 
    s_tm = clock();

    // initialise neighbours
    size_t num_nodes = init_nodes.size(); // initialise number of nodes
    strIntMap init_nbrsize; //
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    for(auto elem: e_id_to_edge){
        auto edge_sz = elem.second.size();
        for(auto v_id: elem.second){
        // initialise number of neighbours and set of neighbours
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp.insert(u);
                        _tmp_sz+=1;
                    }
                }
                init_nbr[v_id] = std::move(_tmp);
                init_nbrsize[v_id] = _tmp_sz;
            }
            else{  // v_id exists in init_nbr map
                auto _tmp = &init_nbr[v_id];
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                init_nbrsize[v_id] = init_nbr[v_id].size();
            }
        }
    }

    // # Initialise buckets
    for (auto node : init_nodes){
        auto len_neighbors = init_nbrsize[node];
        inverse_bucket[node] = len_neighbors;
        if (bucket.find(len_neighbors) == bucket.end()){
            bucket[len_neighbors] = strset();
        }
        bucket[len_neighbors].insert(node);
    }
    clock_t e_tm = clock();
    a.timelogs["init_time"] = double(e_tm - s_tm) / double(CLOCKS_PER_SEC);
    for (size_t k=1; k<= num_nodes; k++){
        while (true){
            if (bucket[k].size()==0)    break;
            auto set_it = bucket[k].begin();  //# get first element in the bucket
            auto v = *set_it;
            bucket[k].erase(set_it);

            a.core[v] = k;
            strvec nbr_v;
            iterate_nbrs(v, nbr_v, inc_dict, e_id_to_edge);
            removeV_transform(v,inc_dict, e_id_to_edge);

            // # enumerating over all neighbors of v
            for (auto u : nbr_v){

                auto len_neighbors_u = get_number_of_nbrs(u, inc_dict, e_id_to_edge);
                if (log)    a.num_nbr_queries += 1;
                // if (nbrquery_stat.find(u) == nbrquery_stat.end()) nbrquery_stat[u] = 0;
                // else    nbrquery_stat[u]+=1;
                auto max_value = std::max(len_neighbors_u, k);

                // # Move u to new location in bucket
                bucket[inverse_bucket[u]].erase(u); // erase u from previous bucket index
                if (bucket.find(max_value) == bucket.end())
                    bucket[max_value] = strset();
                bucket[max_value].insert(u); // insert u to new bucket index
                inverse_bucket[u] = max_value; // update bucket index
            }
        }
    }
    end = clock();
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);

    // if (log){
    //     std::string file = "../output/"+dataset+"_peelnodeQ.csv";
    //     std::stringstream ss2;
    //     std::ofstream out2(file.c_str());
    //     if(out2.fail())
    //     {
    //         out2.close();
    //     }
    //     for(std::string v: init_nodes)
    //     {
    //         ss2<<v<<","<<nbrquery_stat[v]<<"\n";
    //     }
    //     out2 << ss2.str();
    //     out2.close();
    // }
}


// ---------------------------------------------------------------------- Peel ends ---------------------------------------------------------------------------------


// ----------------------------------------------------------------------- E-Peel ------------------------------------------------------------------------------------

void EPeel(std::string dataset, std::map<size_t, strvec > e_id_to_edge, std::map<std::string, std::set<size_t> > inc_dict, strvec init_nodes, Algorithm& a, bool log){
    a.output["algo"] = "E-Peel";
    clock_t start, end;
    // strIntMap nbrquery_stat;
    start = clock();        
    intsStrMap bucket;
    strIntMap inverse_bucket;
    strboolMap setlb;

    clock_t s_tm; 
    s_tm = clock();
    size_t num_nodes = init_nodes.size(); // initialise number of nodes
    strIntMap init_nbrsize; //
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    for(auto elem: e_id_to_edge){
        auto edge_sz = elem.second.size();
        for(auto v_id: elem.second){
        // initialise number of neighbours and set of neighbours
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp.insert(u);
                        _tmp_sz+=1;
                    }
                }
                init_nbr[v_id] = std::move(_tmp);
                init_nbrsize[v_id] = _tmp_sz;
            }
            else{  // v_id exists in init_nbr map
                auto _tmp = &init_nbr[v_id];
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                init_nbrsize[v_id] = init_nbr[v_id].size();
            }
        }
    }
    size_t glb, gub;
    strIntMap lub;
    strIntMap llb;
    glb = std::numeric_limits<size_t>::max();
    gub = std::numeric_limits<size_t>::min();
    // # Computing global lower bounds
    for (auto v :init_nodes){
        auto len_neighbors_v = init_nbrsize[v];
        glb = std::min(glb,len_neighbors_v);
        gub = std::max(gub, len_neighbors_v);
    }
    // Computing Local lower bound 
    for (auto v :init_nodes){
        auto _maxv = std::numeric_limits<size_t>::min();
        for (auto e_id : inc_dict[v]){
            size_t vec_sz = e_id_to_edge[e_id].size();
            _maxv = std::max(_maxv, vec_sz - 1);
        }
        llb[v] = std::max(_maxv, glb);
    }
    auto lb1 = glb;
    auto ub1 = gub;

    // # Initial bucket fill-up
    for (auto node : init_nodes){
        auto lb = llb[node];
        inverse_bucket[node] = lb;
        if (bucket.find(lb) == bucket.end()){
            bucket[lb] = strset();
        }
        bucket[lb].insert(node);
        setlb[node] = true;
    }
    clock_t e_tm = clock();
    a.timelogs["init_time"] = double(e_tm - s_tm) / double(CLOCKS_PER_SEC);

    
    for (size_t k = lb1; k <= ub1; k++){
        while (true){
            if (bucket[k].size()==0)    break;
            auto set_it = bucket[k].begin();  //# get first element in the bucket
            auto v = *set_it;
            bucket[k].erase(set_it);
            if (setlb[v]){
                size_t len_nbr_v = get_number_of_nbrs(v,inc_dict, e_id_to_edge);
                if (log)    a.num_nbr_queries += 1;
                // if (nbrquery_stat.find(v) == nbrquery_stat.end()) nbrquery_stat[v] = 0;
                // else    nbrquery_stat[v]+=1;
                len_nbr_v = std::max(len_nbr_v,k);
                if (bucket.find(len_nbr_v) == bucket.end())
                    bucket[len_nbr_v] = strset();
                bucket[len_nbr_v].insert(v);
                
                // # update new location of u
                inverse_bucket[v] = len_nbr_v;
                setlb[v] = false;
            }
            else{
                a.core[v] = k;
                strvec nbr_v;
                iterate_nbrs(v, nbr_v, inc_dict, e_id_to_edge);
                removeV_transform(v, inc_dict, e_id_to_edge);  //# Store.... + executation time..
    
                for (auto u : nbr_v)
                    if (!setlb[u]){
                        auto len_neighbors_u = get_number_of_nbrs(u, inc_dict, e_id_to_edge);
                        if (log)    a.num_nbr_queries += 1;
                        // if (nbrquery_stat.find(u) == nbrquery_stat.end()) nbrquery_stat[u] = 0;
                        // else    nbrquery_stat[u]+=1;
                        auto  max_value = std::max(len_neighbors_u, k);
                        bucket[inverse_bucket[u]].erase(u);
                        if (bucket.find(max_value) == bucket.end())
                            bucket[max_value] = strset();
                        bucket[max_value].insert(u);
            
    //                     # update new location of u
                        inverse_bucket[u] = max_value;

                    }
            }
        }
    }

    end = clock();
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    
    // if (log){
    //     std::string file = "../output/"+dataset+"_epeellb.csv";
    //     std::stringstream ss;
    //     std::ofstream out(file.c_str());
    //     if(out.fail())
    //     {
    //         std::cout<<"writing failed!\n";
    //         out.close();
    //     }
    //     else{
    //         for(std::string v: init_nodes)
    //         {
    //             size_t max_nbr_Nu = 0;
    //             size_t max_nbr_cu = 0;

    //             size_t num_nbrNu_greater_Nv = 0;
    //             size_t num_nbrNu_greater_cv = 0;
    //             size_t num_nbrNu_greater_lbv = 0;

    //             size_t num_nbrlbu_greater_lbv = 0;

    //             size_t num_nbrcu_greater_cv = 0;
    //             size_t num_nbrcu_noteq_Nv = 0;
    //             size_t num_nbrcu_greater_lbv = 0;
                
    //             for (std::string nbr_v: init_nbr[v]){
    //                     max_nbr_Nu = std::max(init_nbrsize[nbr_v],max_nbr_Nu);
    //                     max_nbr_cu = std::max(a.core[nbr_v],max_nbr_cu);
    //                     if (init_nbrsize[nbr_v] > init_nbrsize[v])  num_nbrNu_greater_Nv += 1;
    //                     if (init_nbrsize[nbr_v] > a.core[v])  num_nbrNu_greater_cv += 1;
    //                     if (init_nbrsize[nbr_v] > llb[v])  num_nbrNu_greater_lbv += 1;

    //                     if (llb[nbr_v] > llb[v])  num_nbrlbu_greater_lbv += 1;

    //                     if (a.core[nbr_v] > a.core[v])  num_nbrcu_greater_cv += 1;
    //                     if (a.core[nbr_v] > llb[v])  num_nbrcu_greater_lbv += 1;
    //                     if (a.core[nbr_v] != init_nbrsize[v])   num_nbrcu_noteq_Nv += 1;
                        
    //             }
    //             // v, lb[v],c[v],|N(v)|, max_{u \in N(v)} |N(u)|, max_{u \in N(v)} c[u]
    //             ss << v << "," << llb[v]<<","<<a.core[v]<<","<<init_nbrsize[v]<<","<<max_nbr_Nu<<","<<max_nbr_cu 
    //             <<","<< num_nbrNu_greater_Nv<<","<< num_nbrNu_greater_cv<<","<< num_nbrNu_greater_lbv
    //             <<","<< num_nbrlbu_greater_lbv<<","<< num_nbrcu_greater_cv <<","<<num_nbrcu_greater_lbv<<","<<num_nbrcu_noteq_Nv<<"\n"; 
    //             // std::cout <<it->first << "," << it->second<<"\n";
    //         }
    //         out << ss.str();
    //         out.close();
    //     }

    //     file = "../output/"+dataset+"_epeelnodeQ.csv";
    //     std::stringstream ss2;
    //     std::ofstream out2(file.c_str());
    //     if(out2.fail())
    //     {
    //         out2.close();
    //     }
    //     for(std::string v: init_nodes)
    //     {
    //         ss2<<v<<","<<nbrquery_stat[v]<<"\n";
    //     }
    //     out2 << ss2.str();
    //     out2.close();
    // }
}

// --------------------------------------------------------------------- E-Peel ends ---------------------------------------------------------------------------------

// --------------------------------------------------------------- Local-Core (Base algorithm) -----------------------------------------------------------------------

bool LCCSAT_check(intvec& inc_edges_u, std::vector<intvec>&edges, size_t u_id, size_t core_u, intvec &hn){

    std::set <size_t> Nplus;
    for (auto e_id :inc_edges_u){
        bool flag = false;
        for(auto v:edges[e_id]){
            if(hn[v]<core_u){
                flag = true;
                break;
            }
        }
        if(!flag){
            for(auto v:edges[e_id]){
                Nplus.insert(v);
            }
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

size_t core_correct(intvec& inc_edges_u, std::vector<intvec>&edges, size_t u_id, size_t core_u, intvec &hn){   // This function is used in optimised local_correct_optIII
        // """ Finds the correct \hat{h} by traversing in descending order from core_u, core_u-1,...,until correct."""
        core_u = core_u - 1; 
        while (LCCSAT_check(inc_edges_u, edges, u_id, core_u, hn) == false) {
            core_u = core_u - 1;
        }

        return core_u;
}

void local_core( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t N = init_nodes.size();
    intvec pcore(N); //
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    strInthashMap node_index; // (use hashtable instead of dictionary => Faster on large |V| datasets.)
    for(size_t i = 0; i<N; i++) node_index[init_nodes[i]] = i; // initialize node_index
    intvec nbrsizes(N); //
    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    std::vector<intvec> nbrs(N);
    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(auto elem: e_id_to_edge){
        auto edge_sz = elem.second.size();
        for(auto v_id: elem.second){
            auto j = node_index[v_id];  // j is of type int
        // initialise number of neighbours and set of neighbours
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp.insert(u);
                        _tmp_sz+=1;
                    }
                }
                init_nbr[v_id] = std::move(_tmp);
                // init_nbrsize[v_id] =  _tmp.size();
                nbrsizes[j] = _tmp_sz;
            }
            else{  // v_id exists in init_nbr map
                auto _tmp = &init_nbr[v_id];
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                // init_nbrsize[v_id] = init_nbr[v_id].size();
                nbrsizes[j] = init_nbr[v_id].size();
            }
        }
    }
    
    end3 = clock();
    time_t start4,end4;
    start4 = clock();
    for (auto node : init_nodes){
        auto _tmp = intvec();
        size_t sz = 0;
        for (auto u: init_nbr[node]){
            _tmp.push_back(node_index[u]);
            sz+=1;
        }
        nbrs[node_index[node]] = std::move(_tmp);
    }
    end4 = clock();
     /* Auxiliary variables for Parallelisation */
    // edge_id is always in [0,M]
    size_t M = 0;
    std::vector< intvec > edges( e_id_to_edge.size() ); // i = edge_id, value = vector of vertices in e[edge_id]
    std::vector<intvec> inc_edges(N); // i=node_id, value = vector of edge ids incident on node_id
    for (auto elem: e_id_to_edge){
        // construct edge
        auto e_id = elem.first;
        auto _tmp = intvec();

        for(auto u: elem.second)  {
            auto i = node_index[u];
            _tmp.push_back(i); 
            inc_edges[i].push_back(e_id);
        }
        edges[M] = std::move(_tmp);
        M+=1;
    }

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrsizes[i]; // initialize pcore
    }
    if (log){
        strstrMap h0;
        for(int i=0;i<N;i++){
            h0[init_nodes[i]] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        // strIntMap h0(pcore);
        a.hnlog.push_back(h0);
    }

    // std::vector<size_t> hn(N);
    intvec hn(N);
    size_t iterations = 0;
    size_t correction_number=0;
    time_t start_main, end_main;
    start_main = clock();
    while (1){
        iterations+=1;
        // std::cout<<"iteration: "<<iteration<<"\n";
        // for(size_t i=0; i<N; i++)
        // {
        //     auto node = hg.init_nodes[i];
        //     std::cout << i<< ": "<< node << "->"<< pcore[i] <<"\n";
        // }
        
        bool flag = true;
        // compute h-index and update core
        for(size_t i = 0; i<N; i++){
            intvec vals(nbrsizes[i]);
            size_t ii=0;
            for (auto u_id : nbrs[i]){
                vals[ii++] = pcore[u_id];
            }
            size_t H_value = hIndex(vals);
            if (H_value < pcore[i]) 
            hn[i] = H_value;
            else hn[i] = pcore[i];
        }
        
        start1 = clock();
        for (size_t i = 0; i<N; i++){
            // auto node = hg.init_nodes[i];
            bool lccsat = LCCSAT_check(inc_edges[i],edges,i,hn[i],hn);
            if (lccsat == false){ 
                flag = false;   //Why this step
                auto hhatn = core_correct(inc_edges[i],edges,i,hn[i],hn);
                pcore[i]  = hhatn;
                // hg.update_min_hindex(node, core[node]);
            }else{
                pcore[i] = hn[i];
            }
        }
        end1 = clock();
        if (log){
            strstrMap h0;
            for(int i=0;i<N;i++){
                h0[init_nodes[i]] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;
    
    }
    end_main = clock();

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        // std::cout << i<< ": "<< node << "->"<< core[node] <<"\n";
    }
    
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);
}

// ------------------------------------------------------------- Local-Core (Base algorithm) ends ---------------------------------------------------------------------

// ----------------------------------------------------------------- Local-Core-OPTI (CSR) ----------------------------------------------------------------------------

bool LCCSAT_check_OPTI(size_t inc_edges_F[], size_t inc_edges_N[], std::vector<intvec>&edges, size_t u_id, size_t core_u, intvec & hn){

    std::set <size_t> Nplus;
    for (size_t i = inc_edges_N[u_id]; i<inc_edges_N[u_id+1]; i++){
        size_t e_id = inc_edges_F[i];
        bool flag = false;
        for(auto v_id: edges[e_id]){
            if(hn[v_id]<core_u){
                flag = true;
                break;
            }
        }
        if (!flag){
            for(auto v_id: edges[e_id]){
                Nplus.insert(v_id);
            }
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

size_t core_correct_OPTI(size_t inc_edges_F[], size_t inc_edges_N[], std::vector<intvec>&edges, size_t u_id, size_t core_u, intvec &hn){   
        // """ Finds the correct \hat{h} by traversing in descending order from core_u, core_u-1,...,until correct."""
        core_u = core_u - 1;
        
        while (LCCSAT_check_OPTI(inc_edges_F,inc_edges_N,edges, u_id, core_u, hn) == false) {
            core_u = core_u - 1;
        }

        return core_u;
}

void local_core_OPTI( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core-OPTI";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t N = init_nodes.size();
    intvec pcore(N); 
    intvec hn(N);

    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    strInthashMap node_index; // (use hashtable instead of dictionary => Faster on large |V| datasets.)

    for(size_t i = 0; i<N; i++) node_index[init_nodes[i]] = i; // initialize node_index
    intvec nbrsizes(N);

    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    std::vector<intvec> nbrs(N);
    
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(auto elem: e_id_to_edge){
        auto edge_sz = elem.second.size();
        sz_inc_edge += edge_sz;
        for(auto v_id: elem.second){
            auto j = node_index[v_id];  
    
        // initialise number of neighbours and set of neighbours
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
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
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                nbrsizes[j] = init_nbr[v_id].size();
            }
        }
    }
    for(auto i : init_nodes){
        sz_init_nbrs += init_nbr[i].size();
    }
    end3 = clock();
    time_t start4,end4;
    start4 = clock();
   
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    nbrs_N[0] = 0;
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));

    int _i = 1, _index=0;
    for (auto node : init_nodes){
        nbrs_N[_i] = nbrs_N[_i-1] + init_nbr[node].size();
        _i++;
        size_t sz = 0;
        for (auto u: init_nbr[node]){
            nbrs_F[_index++] = node_index[u];
            sz+=1;
        }
    }
    end4 = clock();

    size_t M = 0;
    std::vector< intvec > edges( e_id_to_edge.size() ); // i = edge_id, value = vector of vertices in e[edge_id].
    std::vector<intvec> inc_edges(N); // i=node_id, value = vector of edge ids incident on node_id

    for (auto elem: e_id_to_edge){
        // construct edge
        auto e_id = elem.first;
        auto _tmp = intvec();

        for(auto u: elem.second)  {
            auto i = node_index[u];
            _tmp.push_back(i); 
            inc_edges[i].push_back(e_id);
        }
        edges[M] = std::move(_tmp);
        M+=1;
    }
    
    // Calculate csr representation for incident edges
    
    size_t* inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    inc_edges_N[0] = 0;
    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    _i = 1, _index=0;
    for(auto node : init_nodes){
        auto j = node_index[node];
        inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[j].size();
        _i++;
        for(size_t i : inc_edges[j]){
            inc_edges_F[_index++] = i;
        }
    }

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrsizes[i]; // initialize pcore
    }
    
    
    if (log){
        strstrMap h0;
        for(int i=0;i<N;i++){
            h0[init_nodes[i]] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    // std::vector<size_t> hn(N);
    size_t iterations = 0;
    size_t correction_number=0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    start_main = clock();
    while (1){
        iterations+=1;
        
        bool flag = true;
        // compute h-index and update core
        for(size_t i = 0; i<N; i++){
            size_t H_value = hIndex_csr(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if (H_value < pcore[i]) 
                hn[i] = H_value;     
            else
                hn[i] = pcore[i];
        }

        for (size_t i = 0; i<N; i++){
            
            bool lccsat = LCCSAT_check_OPTI(inc_edges_F,inc_edges_N,edges,i,hn[i],hn);
            if (lccsat == false){ 
                flag = false;   //Why this step
                auto hhatn = core_correct_OPTI(inc_edges_F,inc_edges_N,edges,i,hn[i],hn);
                pcore[i]  = hhatn;
                // hg.update_min_hindex(node, core[node]);
            }else{
                pcore[i] = hn[i];
            }
        }
        end1 = clock();
        if (log){
            strstrMap h0;
            for(int i=0;i<N;i++){
                h0[init_nodes[i]] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;
    }
    end_main = clock();
    // std::cout<<"Extra time for calculating N in lccsat_binary "<<t2<<"\n";
    // std::cout<<"Time for h operator "<<hindext<<"\n";
    // std::cout<<"Time for minhindex calculation "<<minht<<"\n";
    // std::cout<<"Extra time for storing hyperedges on the basis of min_hindex "<<t<<"\n";
    // std::cout<<"correction_number "<<correction_number<<"\n";
    // std::cout<<"par_lccsat_count for checking if correction necessary "<<lccsat_count<<"\n";
    // std::cout<<"lccsat_opt_csr_count for core correction "<<lccsat_count2<<"\n";
    // std::cout<<"Time for main loop "<<double(end_main - start_main) / double(CLOCKS_PER_SEC)<<"\n";

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        // std::cout << i<< ": "<< node << "->"<< core[node] <<"\n";
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);
}

// ----------------------------------------------------------------- Local-Core-OPTI (CSR) ends ----------------------------------------------------------------------------

// ----------------------------------------------------------------------- Local-Core-OPTII --------------------------------------------------------------------------------

bool LCCSAT_check_OPTII(size_t inc_edges_F[], size_t inc_edges_N[], std::vector<intvec>&edges, size_t u_id, size_t core_u, intvec & hn){

    std::set <size_t> Nplus;
    for (size_t i = inc_edges_N[u_id]; i<inc_edges_N[u_id+1]; i++){
        size_t e_id = inc_edges_F[i];
        bool flag = false;
        for(auto v_id: edges[e_id]){
            if(hn[v_id]<core_u){
                flag = true;
                break;
            }
        }
        if (!flag){
            for(auto v_id: edges[e_id]){
                Nplus.insert(v_id);
            }
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

size_t core_correct_OPTII(size_t inc_edges_F[], size_t inc_edges_N[], std::vector<intvec>&edges, size_t u_id, size_t core_u, intvec &hn){   
        // """ Finds the correct \hat{h} by traversing in descending order from core_u, core_u-1,...,until correct."""
        core_u = core_u - 1;
        
        while (LCCSAT_check_OPTII(inc_edges_F,inc_edges_N,edges, u_id, core_u, hn) == false) {
            core_u = core_u - 1;
        }

        return core_u;
}

void local_core_OPTII( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core-OPTII";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t N = init_nodes.size();
    intvec pcore(N); 
    // intvec hn(N);

    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    strInthashMap node_index; // (use hashtable instead of dictionary => Faster on large |V| datasets.)

    for(size_t i = 0; i<N; i++) node_index[init_nodes[i]] = i; // initialize node_index
    intvec nbrsizes(N);

    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    std::vector<intvec> nbrs(N);
    
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(auto elem: e_id_to_edge){
        auto edge_sz = elem.second.size();
        sz_inc_edge += edge_sz;
        for(auto v_id: elem.second){
            auto j = node_index[v_id];  
    
        // initialise number of neighbours and set of neighbours
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
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
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                nbrsizes[j] = init_nbr[v_id].size();
            }
        }
    }
    for(auto i : init_nodes){
        sz_init_nbrs += init_nbr[i].size();
    }
    end3 = clock();
    time_t start4,end4;
    start4 = clock();
   
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    nbrs_N[0] = 0;
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));

    int _i = 1, _index=0;
    for (auto node : init_nodes){
        nbrs_N[_i] = nbrs_N[_i-1] + init_nbr[node].size();
        _i++;
        size_t sz = 0;
        for (auto u: init_nbr[node]){
            nbrs_F[_index++] = node_index[u];
            sz+=1;
        }
    }
    end4 = clock();

    size_t M = 0;
    std::vector< intvec > edges( e_id_to_edge.size() ); // i = edge_id, value = vector of vertices in e[edge_id].
    std::vector<intvec> inc_edges(N); // i=node_id, value = vector of edge ids incident on node_id

    for (auto elem: e_id_to_edge){
        // construct edge
        auto e_id = elem.first;
        auto _tmp = intvec();

        for(auto u: elem.second)  {
            auto i = node_index[u];
            _tmp.push_back(i); 
            inc_edges[i].push_back(e_id);
        }
        edges[M] = std::move(_tmp);
        M+=1;
    }
    
    // Calculate csr representation for incident edges
    
    size_t* inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    inc_edges_N[0] = 0;
    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    _i = 1, _index=0;
    for(auto node : init_nodes){
        auto j = node_index[node];
        inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[j].size();
        _i++;
        for(size_t i : inc_edges[j]){
            inc_edges_F[_index++] = i;
        }
    }

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrsizes[i]; // initialize pcore
    }
    
    if (log){
        strstrMap h0;
        for(int i=0;i<N;i++){
            h0[init_nodes[i]] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    // std::vector<size_t> hn(N);
    size_t iterations = 0;
    size_t correction_number=0, check = 0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    start_main = clock();
    while (1){
        iterations+=1;
        
        bool flag = true;
        // compute h-index and update core
        for(size_t i = 0; i<N; i++){
            size_t H_value = hIndex_csr(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if (H_value < pcore[i]) 
                pcore[i] = H_value;     
        }

        for (size_t i = 0; i<N; i++){
            
            check++;
            bool lccsat = LCCSAT_check_OPTII(inc_edges_F,inc_edges_N,edges,i,pcore[i],pcore);
            if (lccsat == false){ 
                correction_number++;
                flag = false;   //Why this step
                auto hhatn = core_correct_OPTII(inc_edges_F,inc_edges_N,edges,i,pcore[i],pcore);
                pcore[i]  = hhatn;
                // hg.update_min_hindex(node, core[node]);
            }
        }
        end1 = clock();
        if (log){
            strstrMap h0;
            for(size_t i=0;i<N;i++){
                h0[init_nodes[i]] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;
    }
    end_main = clock();
    // std::cout<<"Correction number "<<correction_number<<"\n";
    // std::cout<<"Check "<<check<<"\n";
    // std::cout<<"Extra time for calculating N in lccsat_binary "<<t2<<"\n";
    // std::cout<<"Time for h operator "<<hindext<<"\n";
    // std::cout<<"Time for minhindex calculation "<<minht<<"\n";
    // std::cout<<"Extra time for storing hyperedges on the basis of min_hindex "<<t<<"\n";
    // std::cout<<"correction_number "<<correction_number<<"\n";
    // std::cout<<"par_lccsat_count for checking if correction necessary "<<lccsat_count<<"\n";
    // std::cout<<"lccsat_opt_csr_count for core correction "<<lccsat_count2<<"\n";
    // std::cout<<"Time for main loop "<<double(end_main - start_main) / double(CLOCKS_PER_SEC)<<"\n";

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        // std::cout << i<< ": "<< node << "->"<< core[node] <<"\n";
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);
}

// ------------------------------------------------------------------- Local-Core-OPTII ends --------------------------------------------------------------------------------

// ------------------------------------------------------------------- Local-Core-OPTIII -------------------------------------------------------------------------------------

bool LCCSAT_check_OPTIII(size_t inc_edges_F[], size_t inc_edges_N[], std::vector<intvec>&edges, size_t u_id, size_t core_u, intvec & hn){

    std::set <size_t> Nplus;
    for (size_t i = inc_edges_N[u_id]; i<inc_edges_N[u_id+1]; i++){
        size_t e_id = inc_edges_F[i];
        bool flag = false;
        for(auto v_id: edges[e_id]){
            if(hn[v_id]<core_u){
                flag = true;
                break;
            }
        }
        if (!flag){
            for(auto v_id: edges[e_id]){
                Nplus.insert(v_id);
            }
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

size_t core_correct_OPTIII(size_t inc_edges_F[], size_t inc_edges_N[], std::vector<intvec>&edges, size_t u_id, size_t core_u, intvec &hn){   
        // """ Finds the correct \hat{h} by traversing in descending order from core_u, core_u-1,...,until correct."""
        core_u = core_u - 1;
        
        while (LCCSAT_check_OPTIII(inc_edges_F,inc_edges_N,edges, u_id, core_u, hn) == false) {
            core_u = core_u - 1;
        }

        return core_u;
}

void local_core_OPTIII( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core-OPTIII";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t N = init_nodes.size();
    intvec pcore(N); 
    // intvec hn(N);

    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    strInthashMap node_index; // (use hashtable instead of dictionary => Faster on large |V| datasets.)

    for(size_t i = 0; i<N; i++) node_index[init_nodes[i]] = i; // initialize node_index
    intvec nbrsizes(N);

    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    std::vector<intvec> nbrs(N);
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(auto elem: e_id_to_edge){
        auto edge_sz = elem.second.size();
        sz_inc_edge += edge_sz;
        for(auto v_id: elem.second){
            auto j = node_index[v_id];  
            llb[j] = std::max(edge_sz - 1,llb[j]);
        // initialise number of neighbours and set of neighbours
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
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
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                nbrsizes[j] = init_nbr[v_id].size();
            }
        }
    }
    for(auto i : init_nodes){
        sz_init_nbrs += init_nbr[i].size();
    }
    end3 = clock();
    time_t start4,end4;
    start4 = clock();
   
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    nbrs_N[0] = 0;
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));

    size_t glb = std::numeric_limits<size_t>::max();
    size_t gub = std::numeric_limits<size_t>::min();

    int _i = 1, _index=0;
    for (auto node : init_nodes){
        nbrs_N[_i] = nbrs_N[_i-1] + init_nbr[node].size();
        _i++;
        size_t sz = 0;
        for (auto u: init_nbr[node]){
            nbrs_F[_index++] = node_index[u];
            sz+=1;
        }
        glb = std::min(glb, sz);
        gub = std::max(gub, sz);
    }
    end4 = clock();

    size_t M = 0;
    std::vector< intvec > edges( e_id_to_edge.size() ); // i = edge_id, value = vector of vertices in e[edge_id].
    std::vector<intvec> inc_edges(N); // i=node_id, value = vector of edge ids incident on node_id

    for (auto elem: e_id_to_edge){
        // construct edge
        auto e_id = elem.first;
        auto _tmp = intvec();

        for(auto u: elem.second)  {
            auto i = node_index[u];
            _tmp.push_back(i); 
            inc_edges[i].push_back(e_id);
        }
        edges[M] = std::move(_tmp);
        M+=1;
    }
    
    // Calculate csr representation for incident edges
    
    size_t* inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    inc_edges_N[0] = 0;
    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    _i = 1, _index=0;
    for(auto node : init_nodes){
        auto j = node_index[node];
        inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[j].size();
        _i++;
        for(size_t i : inc_edges[j]){
            inc_edges_F[_index++] = i;
        }
    }

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrsizes[i]; // initialize pcore
        llb[i] = std::max(llb[i],glb);
    }
    
    if (log){
        strstrMap h0;
        for(size_t i=0;i<N;i++){
            h0[init_nodes[i]] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    // std::vector<size_t> hn(N);
    size_t iterations = 0;
    size_t correction_number=0, check = 0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    start_main = clock();
    while (1){
        iterations+=1;
        
        bool flag = true;
        // compute h-index and update core
        for(size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            size_t H_value = hIndex_csr(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if (H_value < pcore[i]) 
                pcore[i] = H_value;     
        }

        start1 = clock();
        for (size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            bool lccsat = LCCSAT_check_OPTIII(inc_edges_F,inc_edges_N,edges,i,pcore[i],pcore);
            if (lccsat == false){ 
                start2 = clock();
                correction_number++;
                flag = false;   //Why this step
                auto hhatn = core_correct_OPTIII(inc_edges_F,inc_edges_N,edges,i,pcore[i],pcore);
                pcore[i]  = hhatn;
                end2 = clock();
                a.correction_time+= double(end2-start2)/ double(CLOCKS_PER_SEC);
                // hg.update_min_hindex(node, core[node]);
            }
        }
        end1 = clock();
        a.core_exec_time+= double(end1-start1)/ double(CLOCKS_PER_SEC);
        if (log){
            strstrMap h0;
            for(size_t i=0;i<N;i++){
                h0[init_nodes[i]] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;
    }
    end_main = clock();
    // std::cout<<"Correction number "<<correction_number<<"\n";
    // std::cout<<"Check "<<check<<"\n";
    // std::cout<<"Iterations "<<iterations<<"\n";
    // std::cout<<"Extra time for calculating N in lccsat_binary "<<t2<<"\n";
    // std::cout<<"Time for h operator "<<hindext<<"\n";
    // std::cout<<"Time for minhindex calculation "<<minht<<"\n";
    // std::cout<<"Extra time for storing hyperedges on the basis of min_hindex "<<t<<"\n";
    // std::cout<<"correction_number "<<correction_number<<"\n";
    // std::cout<<"par_lccsat_count for checking if correction necessary "<<lccsat_count<<"\n";
    // std::cout<<"lccsat_opt_csr_count for core correction "<<lccsat_count2<<"\n";
    // std::cout<<"Time for main loop "<<double(end_main - start_main) / double(CLOCKS_PER_SEC)<<"\n";

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        // std::cout << i<< ": "<< node << "->"<< core[node] <<"\n";
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);
}

// ------------------------------------------------------------------- Local-Core-OPTIII ends ---------------------------------------------------------------------------------

// ---------------------------------------------------------------------- Local-Core-OPTIV ------------------------------------------------------------------------------------

bool LCCSAT_check_OPTIV(size_t inc_edges_F[], size_t inc_edges_N[],intvec& min_hindices, std::vector<intvec>&edges, size_t u_id, size_t core_u){

    std::set <size_t> Nplus;
    for (size_t i = inc_edges_N[u_id]; i<inc_edges_N[u_id+1]; i++){
        size_t e_id = inc_edges_F[i];
        if (min_hindices[e_id] >= core_u){
            for(auto v_id: edges[e_id]){
                Nplus.insert(v_id);
            }
            if(Nplus.size()-1>=core_u)
            return true;
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

bool LCCSAT_OPTIV(std::vector<intvec> & min_hindex_to_edge, std::vector<intvec>&edges, size_t core_u, std::set<size_t> & Nplus){

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

size_t core_correct_OPTIV(size_t inc_edges_F[], size_t inc_edges_N[],intvec& min_hindices, std::vector<intvec>&edges, size_t u_id, size_t core_u){   // This function is used in optimised local_correct_optIII
        // """ Finds the correct \hat{h} by traversing in descending order from core_u, core_u-1,...,until correct."""

        // time_t start,end,start1, end1;
        core_u = core_u - 1;

        // start1 = clock();
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
        // end1 = clock();
        std::set<size_t> Nplus;
        
        while (LCCSAT_OPTIV(min_hindex_to_edge, edges, core_u, Nplus) == false) {
            core_u = core_u - 1;
        }

        return core_u;
}

void local_core_OPTIV( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core-OPTIV";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t N = init_nodes.size();
    intvec pcore(N); //
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    strInthashMap node_index; // (use hashtable instead of dictionary => Faster on large |V| datasets.)
    for(size_t i = 0; i<N; i++) node_index[init_nodes[i]] = i; // initialize node_index
    intvec nbrsizes(N); //
    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    std::vector<intvec> nbrs(N);
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(auto elem: e_id_to_edge){
        auto edge_sz = elem.second.size();
        sz_inc_edge += edge_sz;
        for(auto v_id: elem.second){
            auto j = node_index[v_id];  // j is of type int
            llb[j] = std::max(edge_sz - 1,llb[j]);
        // initialise number of neighbours and set of neighbours
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp.insert(u);
                        _tmp_sz+=1;
                    }
                }
                init_nbr[v_id] = std::move(_tmp);
                // init_nbrsize[v_id] =  _tmp.size();
                nbrsizes[j] = _tmp_sz;
                // sz_init_nbrs = _tmp_sz;
            }
            else{  // v_id exists in init_nbr map
                auto _tmp = &init_nbr[v_id];
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                // init_nbrsize[v_id] = init_nbr[v_id].size();
                // sz_init_nbrs -= nbrsizes[j];
                nbrsizes[j] = init_nbr[v_id].size();
                // sz_init_nbrs += nbrsizes[j];
            }
        }
    }
    for(auto i : init_nodes){
        sz_init_nbrs += init_nbr[i].size();
    }
    // std::cout<<"Init nbrs "<<sz_init_nbrs<<" Inc_edges "<<sz_inc_edge<<"\n"; 
    end3 = clock();
    // std::cout<<"Time for init_nbr calculation "<<double(end3 - start3) / double(CLOCKS_PER_SEC)<<"\n";
    time_t start4,end4;
    start4 = clock();

    // size_t nbrs_N[N+1];
    // nbrs_N[0] = 0;
    // size_t nbrs_F[sz_init_nbrs];
   
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    nbrs_N[0] = 0;
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));

    size_t glb = std::numeric_limits<size_t>::max();
    size_t gub = std::numeric_limits<size_t>::min();
    // #pragma omp parallel for default(none) shared(init_nodes,node_index,init_nbr,nbrs)
    int _i = 1, _index=0;
    for (auto node : init_nodes){
        nbrs_N[_i] = nbrs_N[_i-1] + init_nbr[node].size();
        _i++;
        size_t sz = 0;
        for (auto u: init_nbr[node]){
            // _tmp.push_back(node_index[u]);
            nbrs_F[_index++] = node_index[u];
            sz+=1;
        }
        // nbrs[node_index[node]] = std::move(_tmp);
        glb = std::min(glb, sz);
        gub = std::max(gub, sz);
    }
    // nbrs_N[N]++;
    end4 = clock();

    // std::cout<<"Time for nbrs_F (csr-representation) calculation "<<double(end4 - start4) / double(CLOCKS_PER_SEC)<<"\n";
    // std::cout<<"Time for nbrs and init_nbr calculation "<<double(end4 - start3) / double(CLOCKS_PER_SEC)<<"\n";
     /* Auxiliary variables for Parallelisation */
    // edge_id is always in [0,M]
    size_t M = 0;
    std::vector< intvec > edges( e_id_to_edge.size() ); // i = edge_id, value = vector of vertices in e[edge_id]. Try converting this to csr as well
    intvec min_e_hindex( e_id_to_edge.size() ); // i = edge_id, value = minimum of h-index of vertices in e[edge_id]
    std::vector<intvec> inc_edges(N); // i=node_id, value = vector of edge ids incident on node_id

    for (auto elem: e_id_to_edge){
        // construct edge
        auto e_id = elem.first;
        auto _tmp = intvec();
        size_t _min = std::numeric_limits<size_t>::max();

        for(auto u: elem.second)  {
            auto i = node_index[u];
            _tmp.push_back(i); 
            inc_edges[i].push_back(e_id);
            _min = std::min(_min,nbrsizes[i]);
        }
        edges[M] = std::move(_tmp);
        min_e_hindex[M] = _min; // initialize edge h_indices,
        M+=1;
    }
    
    // Calculate csr representation for incident edges
    
    size_t* inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    inc_edges_N[0] = 0;
    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    _i = 1, _index=0;
    for(auto node : init_nodes){
        // std::cout<<_i<<" "<<_index<<"\n";
        auto j = node_index[node];
        inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[j].size();
        _i++;
        for(size_t i : inc_edges[j]){
            inc_edges_F[_index++] = i;
        }
    }

    // initialise core to a upper bound
    // #pragma omp parallel for default(none) shared(nbrsizes,pcore,N,glb,llb)
    for (size_t i = 0; i < N; i++){
        // printf("%d processing: %d\n",omp_get_thread_num(),i);
        pcore[i] = nbrsizes[i]; // initialize pcore
        llb[i] = std::max(llb[i],glb);
    }
    
    if (log){
        strstrMap h0;
        for(size_t i=0;i<N;i++){
            h0[init_nodes[i]] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    // std::vector<size_t> hn(N);
    size_t iterations = 0;
    size_t correction_number=0, check=0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    start_main = clock();
    while (1){
        iterations+=1;
        bool flag = true;
        // compute h-index and update core
        start_h = clock();
        for(size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            size_t H_value = hIndex_csr(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if (H_value < pcore[i]) 
                pcore[i] = H_value;     //pcore[i] is same as hvn here
        }
        end_h = clock();
        hindext += double(end_h - start_h) / double(CLOCKS_PER_SEC);
        // for every edge update is minimum of hindex(constituent vertices)
        start_minh = clock(); 
        for(size_t i = 0; i< M; i++){
            size_t _min = N+1;      //Why M+1 and why calculate 
            for (auto u_id: edges[i])   _min = std::min(_min,pcore[u_id]);
            min_e_hindex[i] = _min;
        }
        end_minh = clock();
        minht += double(end_minh - start_minh) / double(CLOCKS_PER_SEC);
        start1 = clock();
        for (size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            bool lccsat = LCCSAT_check_OPTIV(inc_edges_F,inc_edges_N,min_e_hindex,edges,i,pcore[i]);
            if (lccsat == false){ 
                start2 = clock();  
                correction_number++;
                flag = flag && false;   
                auto hhatn = core_correct_OPTIV(inc_edges_F,inc_edges_N,min_e_hindex,edges,i,pcore[i]);
                pcore[i]  = hhatn;
                for (size_t j = inc_edges_N[i]; j<inc_edges_N[i+1]; j++){
                    size_t e_id = inc_edges_F[j];
                    if (min_e_hindex[e_id] >= pcore[i]){
                        min_e_hindex[e_id] = pcore[i];
                    }
                }
                end2 = clock();
                a.correction_time+= double(end2-start2)/ double(CLOCKS_PER_SEC);
            }
        }
        end1 = clock();
        a.core_exec_time += double(end1 - start1) / double(CLOCKS_PER_SEC);
        end1 = clock();
        if (log){
            strstrMap h0;
            for(size_t i=0;i<N;i++){
                h0[init_nodes[i]] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;

    }
    end_main = clock();
    // std::cout<<"Correction number "<<correction_number<<"\n";
    // std::cout<<"Check "<<check<<"\n";
    // std::cout<<"Iterations "<<iterations<<"\n";
    // std::cout<<"Extra time for calculating N in lccsat_binary "<<t2<<"\n";
    // std::cout<<"Time for h operator "<<hindext<<"\n";
    // std::cout<<"Time for minhindex calculation "<<minht<<"\n";
    // std::cout<<"Extra time for storing hyperedges on the basis of min_hindex "<<t<<"\n";
    // std::cout<<"correction_number "<<correction_number<<"\n";
    // std::cout<<"par_lccsat_count for checking if correction necessary "<<lccsat_count<<"\n";
    // std::cout<<"lccsat_opt_csr_count for core correction "<<lccsat_count2<<"\n";
    // std::cout<<"Time for main loop "<<double(end_main - start_main) / double(CLOCKS_PER_SEC)<<"\n";

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        a.nu_cu += nbrsizes[i] - pcore[i];
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);
}

// ---------------------------------------------------------------------- Local-Core-OPTIV ends --------------------------------------------------------------------------------
