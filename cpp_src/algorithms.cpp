#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include "hypergraph.h"
#include "algorithms.h"
#include "utils.h"

//--------------------------------------------------------------------- Utility functions ------------------------------------------------------------------------------

void Algorithm::printcore(){
    //std::cout << "core: \n";
    for(const auto& elem : core)
    {
    std::cout << elem.first << "->"<<elem.second<<"\n";
    }
}
void Algorithm::writecore(std::string folder){
    //std::cout << "core: \n";
    std::string file = folder + "core_"+output["algo"]+"_"+hg.dataset+".csv";
    std::cout<<"writing to: "<<file<<"\n";
    std::stringstream ss;
    for(auto elem: core)
    {
        // std::cout<<elem.first<<","<<elem.second<<"\n";
        ss << std::to_string(elem.first) << "," << std::to_string(elem.second) << "\n";
    }
    std::ofstream out(file.c_str());
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();
}
void Algorithm::writekdcore(std::string folder){
    std::string file = folder + "core_"+output["algo"]+"_"+hg.dataset+".csv";
    std::cout<<"writing to: "<<file<<"\n";
    std::stringstream ss;
    for(auto elem: core)
    {
        // std::cout<<elem.first<<": "<<elem.second<<","<<secondcore[elem.first]<<"\n";
        ss << std::to_string(elem.first) << "," << std::to_string(elem.second)<<","<<std::to_string(secondcore[elem.first]) << "\n";
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

void iterate_nbrs(size_t v, intvec & nbrs, uintsetvec & inc_dict, intintvec& e_id_to_edge){
    /* Returns the set of neighbours of v.
        implements a traversal from vertex v to each of its neighbours in contrast to set in neighbors(). 
        It also returns an iterator. So it avoids creating the neighborhood list explicitely.
        Overall complexity: O(d(v) * |e_max|), where e_max = largest hyperedge 
    */
    auto incident_edges = inc_dict[v];  //# {O(1)}
    if (incident_edges.size()){
        intboolMap visited_dict;
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

void removeV_transform(size_t v, uintsetvec & inc_dict, intintvec &e_id_to_edge){
    // For every edge e_id incident on v, remove e_id from every vertex u in edge[e_id] distinct from v
    for (auto e_id : inc_dict[v]){
        for(auto u: e_id_to_edge[e_id]){
            if (u==v)   continue;
            inc_dict[u].erase(e_id);
        }
    }
    inc_dict[v] = uintSet();
}

size_t get_number_of_nbrs(size_t v, uintsetvec & inc_dict, intintvec &e_id_to_edge){
    std::unordered_set<int> nbrs;
    auto incident_edges = inc_dict[v];  //# {O(1)}
    if (incident_edges.size()){
        intboolMap visited_dict;
        for (auto e_id : incident_edges){  //# { O(d(v)) }
            for (auto u : e_id_to_edge[e_id]){  //# { O(|e|)}
                if (u != v){
                    nbrs.insert(u);
                }
            }
        }
    }
    return nbrs.size();
}
void degreePeel( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    a.output["algo"] = "naive_deg";
    clock_t start, end;
    start = clock();
    size_t N = init_nodes.size();
    intuSetintMap bucket;
    intvec core(N);
    intvec inverse_bucket(N);
    uintsetvec inc_dict(N, uintSet{});
    intintvec edges( e_id_to_edge.size() ,intvec{});
    for(size_t eid = 0; eid < e_id_to_edge.size(); eid++){
        auto elem = e_id_to_edge[eid];
        auto edge_sz = elem.size();
        for(auto v_id: elem){
            // initialise number of neighbours and set of neighbours
            auto j = node_index[v_id];  
            // std::cout<<j<<" ";
            inc_dict[j].insert(eid);
            edges[eid].push_back(j);
        }
    }
    clock_t e_tm = clock();
    a.output["init_time"] = std::to_string(double(e_tm - start) / double(CLOCKS_PER_SEC));
    start = clock();
    // # Initialise buckets
    size_t max_deg = std::numeric_limits<size_t>::min();
    for (auto node : init_nodes){
        auto i = node_index[node];
        auto degree = inc_dict[i].size();
        max_deg = std::max(max_deg,degree);
        inverse_bucket[i] = degree;
        if (bucket.find(degree) == bucket.end()){
            bucket[degree] = uintSet();
        }
        bucket[degree].insert(i);
    }
    for (size_t k=1; k<= max_deg; k++){
        while (true){
            if (bucket[k].size()==0)    break;
            auto set_it = bucket[k].begin();  //# get first element in the bucket
            auto v = *set_it;
            bucket[k].erase(set_it);
            core[v] = k;
            intvec nbr_v;
            iterate_nbrs(v, nbr_v, inc_dict, edges);
            removeV_transform(v,inc_dict, edges);
            // # enumerating over all neighbors of v
            for (auto u : nbr_v){
                auto degree = inc_dict[u].size();
                auto max_value = std::max(degree, k);

                // # Move u to new location in bucket
                bucket[inverse_bucket[u]].erase(u); // erase u from previous bucket index
                if (bucket.find(max_value) == bucket.end())
                    bucket[max_value] = uintSet();
                bucket[max_value].insert(u); // insert u to new bucket index
                inverse_bucket[u] = max_value; // update bucket index
            }
        }
    }
    end = clock();
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    for(auto node : init_nodes)
    {
        auto i = node_index[node];
        a.core[node] = core[i];
        // std::cout<<a.core[node]<<"-"<<pcore[i]<<"-"<<node<<"-"<<i<<"\n";
    }
}

void Peel( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    a.output["algo"] = "Peel";
    clock_t start, end;
    start = clock();
    size_t N = init_nodes.size();
    intuSetintMap bucket;
    intvec core(N);
    intvec inverse_bucket(N);
    uintsetvec inc_dict(N, uintSet{});
    uintsetvec init_nbr(N, uintSet{});  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    intintvec edges( e_id_to_edge.size() ,intvec{});
    for(size_t eid = 0; eid < e_id_to_edge.size(); eid++){
        auto elem = e_id_to_edge[eid];
        auto edge_sz = elem.size();
        for(auto v_id: elem){
            // initialise number of neighbours and set of neighbours
            auto j = node_index[v_id];  
            // std::cout<<j<<" ";
            inc_dict[j].insert(eid);
            edges[eid].push_back(j);
            auto _tmp = &init_nbr[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
    }
    clock_t e_tm = clock();
    a.output["init_time"] = std::to_string(double(e_tm - start) / double(CLOCKS_PER_SEC));
    start = clock();
    // # Initialise buckets
    for (auto node : init_nodes){
        auto i = node_index[node];
        auto len_neighbors = init_nbr[i].size();
        inverse_bucket[i] = len_neighbors;
        if (bucket.find(len_neighbors) == bucket.end()){
            bucket[len_neighbors] = uintSet();
        }
        bucket[len_neighbors].insert(i);
    }
    for (size_t k=1; k<= N; k++){
        while (true){
            if (bucket[k].size()==0)    break;
            auto set_it = bucket[k].begin();  //# get first element in the bucket
            auto v = *set_it;
            bucket[k].erase(set_it);

            core[v] = k;
            intvec nbr_v;
            iterate_nbrs(v, nbr_v, inc_dict, edges);
            removeV_transform(v,inc_dict, edges);
            // # enumerating over all neighbors of v
            for (auto u : nbr_v){
                auto len_neighbors_u = get_number_of_nbrs(u, inc_dict, edges);
                if (log)    a.num_nbr_queries += 1;
                // if (nbrquery_stat.find(u) == nbrquery_stat.end()) nbrquery_stat[u] = 0;
                // else    nbrquery_stat[u]+=1;
                auto max_value = std::max(len_neighbors_u, k);

                // # Move u to new location in bucket
                bucket[inverse_bucket[u]].erase(u); // erase u from previous bucket index
                if (bucket.find(max_value) == bucket.end())
                    bucket[max_value] = uintSet();
                bucket[max_value].insert(u); // insert u to new bucket index
                inverse_bucket[u] = max_value; // update bucket index
            }
        }
    }
    end = clock();
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    for(auto node : init_nodes)
    {
        auto i = node_index[node];
        a.core[node] = core[i];
        // std::cout<<a.core[node]<<"-"<<pcore[i]<<"-"<<node<<"-"<<i<<"\n";
    }
}

// ----------------------------------------------------------------------- E-Peel ------------------------------------------------------------------------------------
void EPeel( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    a.output["algo"] = "E-Peel";
    clock_t start, end;
    start = clock();
    size_t N = init_nodes.size();
    intuSetintMap bucket;
    intvec core(N);
    std::vector<bool> setlb(N);
    intvec inverse_bucket(N);
    uintsetvec inc_dict(N, uintSet{});
    uintsetvec init_nbr(N, uintSet{});  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    intintvec edges( e_id_to_edge.size() ,intvec{});
    for(size_t eid = 0; eid < e_id_to_edge.size(); eid++){
        auto elem = e_id_to_edge[eid];
        auto edge_sz = elem.size();
        for(auto v_id: elem){
            // initialise number of neighbours and set of neighbours
            auto j = node_index[v_id];  
            // std::cout<<j<<" ";
            inc_dict[j].insert(eid);
            edges[eid].push_back(j);
            auto _tmp = &init_nbr[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
    }
    clock_t e_tm = clock();
    a.output["init_time"] = std::to_string(double(e_tm - start) / double(CLOCKS_PER_SEC));
    start = clock();
    
    size_t glb, gub;
    intvec llb(N);
    glb = std::numeric_limits<size_t>::max();
    gub = std::numeric_limits<size_t>::min();
    // # Computing global lower bounds
    for (size_t i = 0; i<N; i++){
        auto len_neighbors_i = init_nbr[i].size();
        glb = std::min(glb,len_neighbors_i);
        gub = std::max(gub, len_neighbors_i);
    }
        // Computing Local lower bound & bucket init
    for (size_t i=0; i<N; i++){
        auto _maxv = std::numeric_limits<size_t>::min();
        for (auto e_id : inc_dict[i]){
            size_t vec_sz = edges[e_id].size();
            _maxv = std::max(_maxv, vec_sz - 1);
        }
        auto lb = std::max(_maxv, glb);
        llb[i] = lb;
        inverse_bucket[i] = lb;
        if (bucket.find(lb) == bucket.end()){
            bucket[lb] = uintSet();
        }
        bucket[lb].insert(i);
        setlb[i] = true;
    }
    for (size_t k = glb; k <= gub; k++){
        while (true){
            if (bucket[k].size()==0)    break;
            auto set_it = bucket[k].begin();  //# get first element in the bucket
            auto v = *set_it;
            bucket[k].erase(set_it);
            if (setlb[v]){
                size_t len_nbr_v = get_number_of_nbrs(v,inc_dict, edges);
                if (log)    a.num_nbr_queries += 1;
                // if (nbrquery_stat.find(v) == nbrquery_stat.end()) nbrquery_stat[v] = 0;
                // else    nbrquery_stat[v]+=1;
                len_nbr_v = std::max(len_nbr_v,k);
                if (bucket.find(len_nbr_v) == bucket.end())
                    bucket[len_nbr_v] = uintSet();
                bucket[len_nbr_v].insert(v);
                
                // # update new location of u
                inverse_bucket[v] = len_nbr_v;
                setlb[v] = false;
            }
            else{
                core[v] = k;
                intvec nbr_v;
                iterate_nbrs(v, nbr_v, inc_dict, edges);
                removeV_transform(v, inc_dict, edges);  //# Store.... + executation time..
    
                for (auto u : nbr_v)
                    if (!setlb[u]){
                        auto len_neighbors_u = get_number_of_nbrs(u, inc_dict, edges);
                        if (log)    a.num_nbr_queries += 1;
                        // if (nbrquery_stat.find(u) == nbrquery_stat.end()) nbrquery_stat[u] = 0;
                        // else    nbrquery_stat[u]+=1;
                        auto  max_value = std::max(len_neighbors_u, k);
                        bucket[inverse_bucket[u]].erase(u);
                        if (bucket.find(max_value) == bucket.end())
                            bucket[max_value] = uintSet();
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
    for(auto node : init_nodes)
    {
        auto i = node_index[node];
        a.core[node] = core[i];
        // std::cout<<a.core[node]<<"-"<<pcore[i]<<"-"<<node<<"-"<<i<<"\n";
    }
}

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

void local_core( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    size_t N = init_nodes.size();
    intvec pcore(N); //
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)

    start = clock();
    // intvec nbrsizes(N); 
    intintvec edges( e_id_to_edge.size() ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    std::vector<intvec> inc_edges(N, intvec{}); // i=node_id, value = vector of edge ids incident on node_id

    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    // std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    uintsetvec nbrs(N, std::unordered_set<size_t>{});
    // std::vector<intvec> nbrs(N);
    // compute initial neighbors and number of neighbors
    start3 = clock();
    // std::cout<<"sz: "<<e_id_to_edge.size()<<"\n";
    for(size_t eid= 0; eid<e_id_to_edge.size(); eid++){
        // std::cout<<eid<<":\n";
        auto elem = e_id_to_edge[eid];
        // std::cout<<"e: "<<elem.size()<<"\n";
        // auto edge_sz = elem.size();
        for(auto v_id: elem){
            // std::cout<<"vid: "<<v_id<<"-";
            auto j = node_index[v_id];  
            // std::cout<<j<<" ";
            inc_edges[j].push_back(eid);
            edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
    }

    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();
    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
        // std::cout<<i<<":"<<init_nodes[i]<<" "<<pcore[i]<<"\n";
    }
    if (log){
        strstrMap h0;
        for(size_t i=0;i<N;i++){
            h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
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
        // std::cout<<"iteration: "<<iterations<<"\n";
        // for(size_t i=0; i<N; i++)
        // {
        //     auto node = init_nodes[i];
        //     std::cout << i<< ": "<< node << "->"<< pcore[i] <<"\n";
        // }
        
        bool flag = true;
        // compute h-index and update core
        for(size_t i = 0; i<N; i++){
            intvec vals(nbrs[i].size());
            size_t ii=0;
            for (auto u_id : nbrs[i]){
                vals[ii++] = pcore[u_id];
            }
            size_t H_value = hIndex(vals);
            if (H_value < pcore[i]){
                // if (iterations<=2)
                // std::cout<<"updated "<<i<<": "<<pcore[i]<<"=>"<<H_value<<"\n";
                hn[i] = H_value;
            }
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
            for(size_t i=0;i<N;i++){
                h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
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
        // std::cout<<a.core[node]<<"-"<<pcore[i]<<"-"<<node<<"-"<<i<<"\n";
    }
    // std:: cout<<" IT: "<<iterations<<"\n";
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

void local_core_OPTI( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core-OPTI";    
    clock_t start, end, end1;
    // clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices

    size_t N = init_nodes.size();
    intvec pcore(N); //
    intvec hn(N);
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    intintvec edges( e_id_to_edge.size() ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    intintvec inc_edges(N, intvec{}); // i=node_id, value = vector of edge ids incident on node_id

    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    // std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    uintsetvec nbrs(N, std::unordered_set<size_t>{});
    // compute initial neighbors and number of neighbors
    for(size_t eid= 0; eid<e_id_to_edge.size(); eid++){
        auto elem = e_id_to_edge[eid];
        sz_inc_edge += elem.size();
        for(auto v_id: elem){
            auto j = node_index[v_id];
            inc_edges[j].push_back(eid);
            edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
        // std::cout<<"\n";
    }
    for(size_t _i = 0; _i< N; _i ++){
        sz_init_nbrs += nbrs[_i].size();
        // std::cout<<init_nodes[_i]<<" nbrs: ";
        // for(auto u:nbrs[_i])    std::cout<<init_nodes[u]<<",";
        // std::cout<<"\n";
    }
    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    size_t *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
    inc_edges_N[0] = 0;
    nbrs_N[0] = 0;
    for (int _i = 1; _i<= N; _i ++){
        nbrs_N[_i] = nbrs_N[_i-1] + nbrs[_i-1].size();
		inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[_i-1].size();
    }
    
    // Calculate csr representation for incident edges
    for (int _i = 1; _i<= N; _i ++){
		auto _index = nbrs_N[_i-1];
		for(auto u: nbrs[_i-1]){
			nbrs_F[_index++] = u;
		}
		_index = inc_edges_N[_i-1];
		for(auto eid : inc_edges[_i-1])
			inc_edges_F[_index++] = eid;
	}
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();
    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
    }
    
    if (log){
        strstrMap h0;
        for(int i=0;i<N;i++){
            h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    size_t iterations = 0;
    size_t correction_number=0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    // start_main = clock();
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
            }else{
                pcore[i] = hn[i];
            }
        }
        end1 = clock();
        if (log){
            strstrMap h0;
            for(int i=0;i<N;i++){
                h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;
    }
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

void local_core_OPTII(std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core-OPTII";
    clock_t start, end, end1;
    // clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices

    size_t N = init_nodes.size();
    intvec pcore(N); //
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    intintvec edges( e_id_to_edge.size() ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    intintvec inc_edges(N, intvec{}); // i=node_id, value = vector of edge ids incident on node_id

    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    // std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    uintsetvec nbrs(N, std::unordered_set<size_t>{});
    // compute initial neighbors and number of neighbors
    for(size_t eid= 0; eid<e_id_to_edge.size(); eid++){
        auto elem = e_id_to_edge[eid];
        sz_inc_edge += elem.size();
        for(auto v_id: elem){
            auto j = node_index[v_id];
            inc_edges[j].push_back(eid);
            edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
        // std::cout<<"\n";
    }
    for(size_t _i = 0; _i< N; _i ++){
        sz_init_nbrs += nbrs[_i].size();
        // std::cout<<init_nodes[_i]<<" nbrs: ";
        // for(auto u:nbrs[_i])    std::cout<<init_nodes[u]<<",";
        // std::cout<<"\n";
    }
    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    size_t *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
    inc_edges_N[0] = 0;
    nbrs_N[0] = 0;
    for (int _i = 1; _i<= N; _i ++){
        nbrs_N[_i] = nbrs_N[_i-1] + nbrs[_i-1].size();
		inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[_i-1].size();
    }
    
    // Calculate csr representation for incident edges
    for (int _i = 1; _i<= N; _i ++){
		auto _index = nbrs_N[_i-1];
		for(auto u: nbrs[_i-1]){
			nbrs_F[_index++] = u;
		}
		_index = inc_edges_N[_i-1];
		for(auto eid : inc_edges[_i-1])
			inc_edges_F[_index++] = eid;
	}
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();
    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
    }
    
    if (log){
        strstrMap h0;
        for(int i=0;i<N;i++){
            h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    size_t iterations = 0;
    size_t correction_number=0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    // start_main = clock();
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
            bool lccsat = LCCSAT_check_OPTII(inc_edges_F,inc_edges_N,edges,i,pcore[i],pcore);
            if (lccsat == false){ 
                flag = false;   //Why this step
                auto hhatn = core_correct_OPTII(inc_edges_F,inc_edges_N,edges,i,pcore[i],pcore);
                pcore[i]  = hhatn;
            }
        }
        end1 = clock();
        if (log){
            strstrMap h0;
            for(int i=0;i<N;i++){
                h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;
    }
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

void local_core_OPTIII( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log){   
    a.output["algo"] = "Local-core-OPTIII";
    clock_t start, end, end1;
    // clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    size_t N = init_nodes.size();
    intvec pcore(N); //
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t glb = std::numeric_limits<size_t>::max();

    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    intintvec edges( e_id_to_edge.size() ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    intintvec inc_edges(N, intvec{}); // i=node_id, value = vector of edge ids incident on node_id

    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    // std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    uintsetvec nbrs(N, std::unordered_set<size_t>{});
    // compute initial neighbors and number of neighbors
    for(size_t eid= 0; eid<e_id_to_edge.size(); eid++){
        auto elem = e_id_to_edge[eid];
        sz_inc_edge += elem.size();
        for(auto v_id: elem){
            auto j = node_index[v_id];
            inc_edges[j].push_back(eid);
            edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
        // std::cout<<"\n";
    }
    for(size_t _i = 0; _i< N; _i ++){
        sz_init_nbrs += nbrs[_i].size();
        // std::cout<<init_nodes[_i]<<" nbrs: ";
        // for(auto u:nbrs[_i])    std::cout<<init_nodes[u]<<",";
        // std::cout<<"\n";
    }
    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    size_t *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
    inc_edges_N[0] = 0;
    nbrs_N[0] = 0;
    for (int _i = 1; _i<= N; _i ++){
        auto nbr_i = nbrs[_i-1].size();
        nbrs_N[_i] = nbrs_N[_i-1] + nbr_i;
        glb = std::min(glb, nbr_i);
		inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[_i-1].size();
    }
    
    // Calculate csr representation for incident edges
    for (int _i = 1; _i<= N; _i ++){
		auto _index = nbrs_N[_i-1];
		for(auto u: nbrs[_i-1]){
			nbrs_F[_index++] = u;
		}
		_index = inc_edges_N[_i-1];
		for(auto eid : inc_edges[_i-1])
			inc_edges_F[_index++] = eid;
	}
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();
    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
        llb[i] = std::max(llb[i],glb);
    }
    
    if (log){
        strstrMap h0;
        for(int i=0;i<N;i++){
            h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    size_t iterations = 0;
    size_t correction_number=0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    // start_main = clock();
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

        for (size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            bool lccsat = LCCSAT_check_OPTIII(inc_edges_F,inc_edges_N,edges,i,pcore[i],pcore);
            if (lccsat == false){ 
                flag = false;   //Why this step
                auto hhatn = core_correct_OPTIII(inc_edges_F,inc_edges_N,edges,i,pcore[i],pcore);
                pcore[i]  = hhatn;
            }
        }
        end1 = clock();
        if (log){
            strstrMap h0;
            for(int i=0;i<N;i++){
                h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;
    }
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

void local_core_OPTIV( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    a.output["algo"] = "Local-core-OPTIV";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    size_t N = init_nodes.size();
    size_t M = e_id_to_edge.size();
    intvec pcore(N); //
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t glb = std::numeric_limits<size_t>::max();
    intintvec edges( M ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    intvec min_e_hindex(M);
    intintvec inc_edges(N, intvec{}); // i=node_id, value = vector of edge ids incident on node_id
    uintsetvec nbrs(N, std::unordered_set<size_t>{});

    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(size_t eid= 0; eid<M; eid++){
        auto elem = e_id_to_edge[eid];
        sz_inc_edge += elem.size();
        for(auto v_id: elem){
            auto j = node_index[v_id];
            inc_edges[j].push_back(eid);
            edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
        // std::cout<<"\n";
    }
    for(size_t _i = 0; _i< N; _i ++){
        sz_init_nbrs += nbrs[_i].size();
    }

    // std::cout<<"Init nbrs "<<sz_init_nbrs<<" Inc_edges "<<sz_inc_edge<<"\n"; 
    end3 = clock();
    // std::cout<<"Time for init_nbr calculation "<<double(end3 - start3) / double(CLOCKS_PER_SEC)<<"\n";
    time_t start4,end4;
    start4 = clock();

    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    size_t *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
    inc_edges_N[0] = 0;
    nbrs_N[0] = 0;
    for (int _i = 1; _i<= N; _i ++){
        auto nbr_i = nbrs[_i-1].size();
        nbrs_N[_i] = nbrs_N[_i-1] + nbr_i;
        glb = std::min(glb, nbr_i);
		inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[_i-1].size();
    }
    
    // Calculate csr representation for incident edges
    for (int _i = 1; _i<= N; _i ++){
		auto _index = nbrs_N[_i-1];
		for(auto u: nbrs[_i-1]){
			nbrs_F[_index++] = u;
		}
		_index = inc_edges_N[_i-1];
		for(auto eid : inc_edges[_i-1])
			inc_edges_F[_index++] = eid;
	}
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
        llb[i] = std::max(llb[i],glb);
    }
    
    if (log){
        strstrMap h0;
        for(int i=0;i<N;i++){
            h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
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
                h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
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
        a.nu_cu += nbrs[i].size() - pcore[i];
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);

}

// ---------------------------------------------------------------------- Local-Core-OPTIV ends --------------------------------------------------------------------------------

// ------------------------------------------------------- Local-core Clique graph -------------------------
// void local_core_clique( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a){   
void local_core_clique( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    std::cout<<"clique-core decomposition\n";
    log = false;
     a.output["algo"] = "clique";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    // size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    size_t N = init_nodes.size();
    size_t M = e_id_to_edge.size();
    intvec pcore(N); //
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t glb = std::numeric_limits<size_t>::max();
    // intintvec edges( M ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    // intvec min_e_hindex(M);
    // intintvec inc_edges(N, intvec{}); // i=node_id, value = vector of edge ids incident on node_id
    uintsetvec nbrs(N, std::unordered_set<size_t>{});

    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(size_t eid= 0; eid<M; eid++){
        auto elem = e_id_to_edge[eid];
        // sz_inc_edge += elem.size();
        for(auto v_id: elem){
            auto j = node_index[v_id];
            // inc_edges[j].push_back(eid);
            // edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
        // std::cout<<"\n";
    }
    for(size_t _i = 0; _i< N; _i ++){
        sz_init_nbrs += nbrs[_i].size();
    }

    // std::cout<<"Init nbrs "<<sz_init_nbrs<<" Inc_edges "<<sz_inc_edge<<"\n"; 
    end3 = clock();
    // std::cout<<"Time for init_nbr calculation "<<double(end3 - start3) / double(CLOCKS_PER_SEC)<<"\n";
    time_t start4,end4;
    start4 = clock();

    // size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    // size_t *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
    // inc_edges_N[0] = 0;
    nbrs_N[0] = 0;
    for (int _i = 1; _i<= N; _i ++){
        auto nbr_i = nbrs[_i-1].size();
        nbrs_N[_i] = nbrs_N[_i-1] + nbr_i;
        glb = std::min(glb, nbr_i);
		// inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[_i-1].size();
    }
    
    // Calculate csr representation for incident edges
    for (int _i = 1; _i<= N; _i ++){
		auto _index = nbrs_N[_i-1];
		for(auto u: nbrs[_i-1]){
			nbrs_F[_index++] = u;
		}
		// _index = inc_edges_N[_i-1];
		// for(auto eid : inc_edges[_i-1])
		// 	inc_edges_F[_index++] = eid;
	}
    end = clock();
    a.output["init_time"] = std::to_string((double(end - start) / double(CLOCKS_PER_SEC)));
    start = clock();

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
        llb[i] = std::max(llb[i],glb);
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
        for(size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            size_t H_value = hIndex_csr(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if (H_value != pcore[i])    flag = false; 
            pcore[i] = H_value;     //pcore[i] is same as hvn here
        }
        end1 = clock();
        if (flag)
            break;

    }
    end_main = clock();
    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        a.nu_cu += nbrs[i].size() - pcore[i];
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);

}

void print_bucket(intuSetintMap& degbucket, intvec& init_nodes){
    for(auto pr: degbucket) {
        std::cout<<pr.first<<": ";
        for(auto u: pr.second) std::cout<< init_nodes[u]<<", "; 
        std::cout<<"\n";
    }
}
void kdCorehybrid(std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    a.output["algo"] = "kdcore";
    clock_t start, end;
    double init_tm=0, actual_tm=0;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    size_t N = init_nodes.size();
    size_t M = e_id_to_edge.size();
    intvec pcore(N); //
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t glb = std::numeric_limits<size_t>::max();
    intintvec edges( M ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    intvec min_e_hindex(M);
    uintsetvec inc_edges(N, uintSet{});
    uintsetvec nbrs(N, std::unordered_set<size_t>{});
    std::vector<std::vector<intpair>> score_mult(N,std::vector<intpair>{});
    // compute initial neighbors and number of neighbors
    for(size_t eid= 0; eid<M; eid++){
        auto elem = e_id_to_edge[eid];
        sz_inc_edge += elem.size();
        for(auto v_id: elem){
            auto j = node_index[v_id];
            inc_edges[j].insert(eid);
            edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
    }
    for(size_t _i = 0; _i< N; _i ++){
        sz_init_nbrs += nbrs[_i].size();
    }

    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    size_t *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
    inc_edges_N[0] = 0;
    nbrs_N[0] = 0;
    for (int _i = 1; _i<= N; _i ++){
        auto nbr_i = nbrs[_i-1].size();
        nbrs_N[_i] = nbrs_N[_i-1] + nbr_i;
        glb = std::min(glb, nbr_i);
		inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[_i-1].size();
    }
    
    // Calculate csr representation for incident edges
    for (int _i = 1; _i<= N; _i ++){
		auto _index = nbrs_N[_i-1];
		for(auto u: nbrs[_i-1]){
			nbrs_F[_index++] = u;
		}
		_index = inc_edges_N[_i-1];
		for(auto eid : inc_edges[_i-1])
			inc_edges_F[_index++] = eid;
	}
    init_tm = (double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
        llb[i] = std::max(llb[i],glb);
    }

    // std::vector<size_t> hn(N);
    size_t iterations = 0;
    while (1){
        iterations+=1;
        bool flag = true;
        // compute h-index and update core
        for(size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            size_t H_value = hIndex_csr(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if (H_value < pcore[i]) 
                pcore[i] = H_value;     //pcore[i] is same as hvn here
        }
        for(size_t i = 0; i< M; i++){
            size_t _min = N+1;      //Why M+1 and why calculate 
            for (auto u_id: edges[i])   _min = std::min(_min,pcore[u_id]);
            min_e_hindex[i] = _min;
        }
        for (size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            bool lccsat = LCCSAT_check_OPTIV(inc_edges_F,inc_edges_N,min_e_hindex,edges,i,pcore[i]);
            if (lccsat == false){ 
                flag = flag && false;   
                auto hhatn = core_correct_OPTIV(inc_edges_F,inc_edges_N,min_e_hindex,edges,i,pcore[i]);
                pcore[i]  = hhatn;
                for (size_t j = inc_edges_N[i]; j<inc_edges_N[i+1]; j++){
                    size_t e_id = inc_edges_F[j];
                    if (min_e_hindex[e_id] >= pcore[i]){
                        min_e_hindex[e_id] = pcore[i];
                    }
                }
            }
        }
        if (flag)
            break;

    }
    // Peeling iteration to find secondary core.
    // intuSetintMap nbrbucket;
    intvec inverse_bucket(N);
	// initialize every nodes initial bucket to the primary core-number.
	size_t min_cv = N+1;
    size_t max_cv = 0;
    for (size_t i = 0; i<N; i++){
        auto cv = pcore[i];
        // auto v = i;
		// if (nbrbucket.find(cv) == nbrbucket.end())
		// 	nbrbucket[cv] = uintSet({v});
		// else
		// 	nbrbucket[cv].insert(v);
        max_cv = std::max(max_cv, cv);
        min_cv = std::min(min_cv, cv);
    }
    // std::cout<<min_cv<<":"<<max_cv<<"\n";
	for (size_t pk = min_cv; pk<= max_cv; pk++){
        intvec score(N,-1); //
        if(log) std::cout<<"pcore="<<pk<<"\n";
		// deg bucket init 
		intuSetintMap degbucket;
		size_t max_deg = 0;
        std::vector<bool>stop(N,false);
		for (size_t u = 0; u<N; u++){
            if(pcore[u]>=pk){
                auto d = inc_edges[u].size();
                if (degbucket.find(d) == degbucket.end()) degbucket[d] = uintSet();
                degbucket[d].insert(u);
                inverse_bucket[u] = d;
                max_deg = std::max(d,max_deg);
            }
		}
        if(log) {
            std::cout<<"max_deg: "<<max_deg<<"\n";
            std::cout<<"init degbucket: \n"; 
            print_bucket(degbucket,init_nodes);
            std::cout<<"inc_edges: \n";
            for(size_t i =0; i<N; i++){
                std::cout<<init_nodes[i]<<": ";
                for(auto u: inc_edges[i])  std::cout<<init_nodes[u]<<", ";
                std::cout<<"\n";
            }
        }
		// bool stop = false;
		// size_t maximal_dk = 1;
		for(size_t dk = 1; dk<= max_deg; dk++){
            if(log) std::cout<<"dk = "<<dk<<"\n";
            if(degbucket.find(dk)==degbucket.end()) continue;
			while (degbucket[dk].size()!=0){
				// Pop v from degbucket[dk];
                auto set_it = degbucket[dk].begin();  //# get first element in the bucket
                auto v = *set_it;
                degbucket[dk].erase(set_it);
                if (log){
                    std::cout<<"pop: "<<v<<"/"<<init_nodes[v]<<"\n";
                }
				score[v] = dk; // assign secondary core-num to v
                stop[v] = true;
                intvec nbrs_v;
                iterate_nbrs(v, nbrs_v, inc_edges, edges);
                if(log){
                    std::cout<<"iterate_nbrs: \n";
                    for(auto u: nbrs_v) std::cout<<u<<"/"<<init_nodes[u]<<",";
                    std::cout<<"\n";
                }
                removeV_transform(v,inc_edges, edges);
                if(log){
                    std::cout<<"inc_edge after removal: \n";
                    for(size_t i =0; i<N; i++){
                        std::cout<<init_nodes[i]<<": ";
                        for(auto u: inc_edges[i])  std::cout<<init_nodes[u]<<", ";
                        std::cout<<"\n";
                    }
                    std::cout<<"done\n";
                }
				for(auto u: nbrs_v){
                    if(log) std::cout<<" -- "<<u<<"/"<<init_nodes[u]<<"\n"; 
                    if(stop[u]) continue;
					// if |N(u)| in residual hyp < primary core , stop 
					if (get_number_of_nbrs(u, inc_edges, edges)< pk){
						// stop  peeling v caused nbr u's |N(u)| in the current subhyp. < pk
                        if(log) std::cout<<"stop\n"; 
                        if(log){
                            std::cout<<"current deg: \n";
                            for(size_t i = 0; i<N; i++){
                                std::cout<<init_nodes[i]<<": "<<inc_edges[i].size()<<"\n";
                            }
                        }
                        degbucket[inverse_bucket[u]].erase(u); // erase u from previous bucket index
                        // if (degbucket.find(dk) == degbucket.end()) degbucket[dk] = uintSet();
                        degbucket[dk].insert(u);

                        if (log) {std::cout<<"bucket: \n"; print_bucket(degbucket,init_nodes);}
                        if(log) std::cout<<"batch delete\n";
                        /* We peel remaining nodes with pcore[u] == pk one by one without 
                        doing expensive nbr traversal for efficiency. 
                        We could have taken induced subhypergrpah {u: pcore[u]>=pk} but that would 
                        require constructing sub-hyp. from scratch which is again more expensive than just
                        peeling the remainder nodes.
                        */
                        // for(size_t ddk = dk; ddk<=max_deg; ddk++){
                        //     if (log) std::cout<<"ddk: "<<ddk<<"\n";
                        //     if(degbucket.find(ddk)!=degbucket.end()){
                        //         while(degbucket[ddk].size()){
                        //             auto set_it = degbucket[ddk].begin();  //# get first element in the bucket
                        //             auto u = *set_it;
                        //             if (log) std::cout<<init_nodes[u]<<",";
                        //             degbucket[ddk].erase(set_it);
                        //             if(pcore[u]==pk)
                        //                 removeV_transform(u,inc_edges, edges);
                        //             score[u] = dk;
                        //         }
                        //     }
                        //     if(log) std::cout<<"\n";
                        // }
                        // break;
					}
					else{ // else, update index in degree bucket for u \in N(v) 
                        // only update bucket position of nodes in nbr pk-core.
                        // pk+1, and higher core-nodes will be processed in later time.
                        // if (nbrbucket[pk].find(u) != nbrbucket[pk].end()){
                            auto d = inc_edges[u].size();
                            d = std::max(d,dk);
                            degbucket[inverse_bucket[u]].erase(u); // erase u from previous bucket index
                            if (degbucket.find(d) == degbucket.end()) degbucket[d] = uintSet();
                            degbucket[d].insert(u);
                            if (log){std::cout<< "bucket update: \n";    print_bucket(degbucket,init_nodes);}
                            inverse_bucket[u] = d;
                        // }
					}
				}
                if (log) std::cout<<"done traversing nbrs\n";
			}
		}
        for(size_t i=0; i<N; i++){
            if(score[i]!=-1)
                score_mult[i].push_back(std::make_pair(pcore[i],score[i]));
        }
    }
    a.exec_time = double(clock() - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);
    a.output["init_time"] = std::to_string(init_tm);
    // for(size_t i=0; i<N; i++){
    //     auto node = init_nodes[i];
    //     a.core[node] = pcore[i];
    //     a.secondcore[node] = score[i];
    //     // std::cout << i<< ": "<< node << "->"<< pcore[i] <<" , "<<score[i]<<"\n";
    // }
    for(size_t i = 0; i<N; i++){
        auto node = init_nodes[i];
        std::cout<<node<<": ";
        auto cores = score_mult[i].size();
        for(size_t j = 0; j<cores; j++){
            std::cout<<score_mult[i][j].first<<","<<score_mult[i][j].second<<"\t";
        }
        std::cout<<"\n";
    }
}

// ------------------------------------------------------- Local-core dist2 Bipartite graph -------------------------

void getBipartite(intintvec &e_id_to_edge, std::unordered_map<size_t,intvec> &v_to_eid, std::unordered_map<size_t,intvec> &eid_to_v){

    size_t index = 0;
    for(auto edge:e_id_to_edge){
        
        for(auto v:edge){
            v_to_eid[v].push_back(index);
            eid_to_v[index].push_back(v);
        }
        index++;
    }

}

std::set<size_t> getDist2nbr(std::unordered_map<size_t,intvec> &v_to_eid, std::unordered_map<size_t,intvec> &eid_to_v, size_t u){

    std::set<size_t> s;
    for(auto eid:v_to_eid[u]){
        for(auto v:eid_to_v[eid]){
            if(v!=u)
            s.insert(v);
        }
    }
    return s;

}

std::set<int> getDist2nbredge(std::unordered_map<size_t,intvec> &v_to_eid, std::unordered_map<size_t,intvec> &eid_to_v, size_t u){

    std::set<int> s;
    for(auto v:eid_to_v[u]){
        for(auto e:v_to_eid[v]){
            if(e!=u)
            s.insert(e);
        }
    }
    return s;

}

void local_core_bipartite(std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    a.output["algo"] = "bipartite";
    time_t start = clock();
    std::unordered_map<size_t,size_t> core;
    int num_edge = e_id_to_edge.size();
    for(auto i:node_index){
        node_index[i.first] += num_edge;
    }
    std::unordered_map<size_t,intvec> v_to_eid;
    std::unordered_map<size_t,intvec> eid_to_v;
    getBipartite(e_id_to_edge,v_to_eid,eid_to_v);
    intvec pcore(eid_to_v.size()+v_to_eid.size());

    //Initialise pcore
    for(size_t v:init_nodes){
        pcore[node_index[v]] = getDist2nbr(v_to_eid,eid_to_v,v).size() + v_to_eid[v].size();
        // std::cout<<v<<" "<<pcore[node_index[v]]<<"\n";
    }
    // std::cout<<"edges\n";
    for(int e_id=0;e_id<e_id_to_edge.size();e_id++){
        pcore[e_id] = getDist2nbredge(v_to_eid,eid_to_v,e_id).size() + eid_to_v[e_id].size();
        // std::cout<<e.first<<" "<<pcore[e.first]<<"\n";
    }
    // end = clock();
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();
    int iterations = 0;
    while(1){
        // std::cout<<iterations<<"\n";
        iterations++;
        // for(int i=0;i<pcore.size();i++)
        // std::cout<<pcore[i]<<" ";
        // std::cout<<"\n";
        // assert(eid_to_v.size()==e_id_to_edge.size());
        // assert(v_to_eid.size()==init_nodes.size());
        intvec H(eid_to_v.size()+v_to_eid.size());
        for(auto v:init_nodes){

            //Given H(n-1), calculate A(n-1). 
            // attainability is per start node and end node with distance <= 2 from start node. 

            std::unordered_map<size_t,size_t> attainability;
            for(auto e:v_to_eid[v]){
                attainability[e] = pcore[e];
                for(auto u:eid_to_v[e]){
                    if(u!=v){
                        if(attainability.find(node_index[u])==attainability.end()){
                            attainability[node_index[u]] = std::min(pcore[e],pcore[node_index[u]]);
                        }else{
                            size_t curr_attain = std::min(pcore[e],pcore[node_index[u]]);
                            attainability[node_index[u]] = std::max(attainability[node_index[u]],curr_attain);
                        }
                    }
                }
            }

            // Use A(n-1) to calculate H(n).
            intvec vals;
            for(auto val : attainability)
            vals.push_back(val.second);
            H[node_index[v]] = hIndex(vals);

        }

        for(size_t e=0;e<e_id_to_edge.size();e++){

            //Given H(n-1), calculate A(n-1). 
            // attainability is per start node and end node with distance <= 2 from start node. 

            std::unordered_map<size_t,size_t> attainability;
            for(auto v:eid_to_v[e]){
                attainability[node_index[v]] = pcore[node_index[v]];
                for(auto eid : v_to_eid[v]){
                    if(e!=eid){
                        if(attainability.find(eid)==attainability.end()){
                            attainability[eid] = std::min(pcore[eid],pcore[node_index[v]]);
                        }else{
                            size_t curr_attain = std::min(pcore[eid],pcore[node_index[v]]);
                            attainability[eid] = std::max(attainability[eid],curr_attain);
                        }
                    }
                }
            }

            // Use A(n-1) to calculate H(n).
            intvec vals;
            for(auto val : attainability)
            vals.push_back(val.second);
            H[e] = hIndex(vals);

        }

        bool stop = true;
        int count = 0;
        for(int i=0;i<H.size();i++){
            if(pcore[i]!=H[i]){
                count++;
                stop = false;
            }
            pcore[i] = H[i];
        }
        // std::cout<<count<<"\n";

        if(stop)
        break;

    }
    a.exec_time = double(clock() - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);

    for(int i:init_nodes){
        a.core[i] = pcore[node_index[i]];
    }
}