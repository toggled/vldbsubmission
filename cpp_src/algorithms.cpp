#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include <algorithm>
#include <unordered_set>
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
    std::cout << "core: \n";
    std::string file = "../output/core_"+output["algo"]+"_"+hg.dataset+".csv";
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

void iterate_nbrs(size_t v, intvec & nbrs, uintsetvec & inc_dict, intintvec e_id_to_edge){
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
