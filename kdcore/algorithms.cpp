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
    std::string file = "../output/"+output["algo"]+"_"+hg.dataset+".csv";
    std::stringstream ss;
    for(auto node: hg.init_nodes){
         ss << node <<","<<std::to_string(core[node])<<","<<std::to_string(secondcore[node])<<"\n";
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
    std::string file = "../output/kdresults.csv";
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

void  Algorithm::writelog(){
    std::string file = "../output/kdlog_"+output["algo"]+"_"+hg.dataset+".csv";
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

void local_kdcore( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a, bool log){
    a.output["algo"] = "kdcore";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t N = init_nodes.size();
    intvec pcore(N); //
    intvec score(N); //
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
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrsizes[i]; // initialize pcore
        score[i] =  inc_edges[i].size();
    }
    if (log){
        strstrMap h0;
        // strstrMap g0;
        for(int i=0;i<N;i++){
            h0[init_nodes[i]] = std::to_string(pcore[i]);
            // g0[init_nodes[i]] = std::to_string(score[i]);
        }
        h0["Time"] = "0";
        // strIntMap h0(pcore);
        a.hnlog.push_back(h0);
        // a.gnlog.push_back(g0);
    }
    intvec hn(N);
    intvec gn(N);
    size_t iterations = 0;
    size_t correction_number=0;
    time_t start_main, end_main;
    start_main = clock();
    while (1){
        iterations+=1;
        
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
            bool lccsat = LCCSAT_check(inc_edges[i],edges,i,hn[i],hn);
            if (lccsat == false){ 
                flag = false;   //Why this step
                auto hhatn = core_correct(inc_edges[i],edges,i,hn[i],hn);
                pcore[i]  = hhatn;
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
    a.output["total iteration(p)"] = std::to_string(iterations);

    // Compute secondary core-numbers
    iterations = 0;
    intIntMap eid_minh;
    intIntMap eid_ming;
    if (log){
        std::cout<<"primary core: \n";
        for(size_t i=0; i<N; i++)
        {
            std::cout<<init_nodes[i]<<": "<<pcore[i]<<"\n";
        }
    }

    // Secondary-core computation loop
    while (1){
        iterations+=1;
        bool flag = true;
        for (auto elem: e_id_to_edge){
            auto e_id = elem.first;
            auto _tmp = intvec();
            size_t _mins = std::numeric_limits<size_t>::max();
            size_t _minp = std::numeric_limits<size_t>::max();
            bool same_pcore = true;
            for(auto u: elem.second)  {
                auto i = node_index[u];
                _mins = std::min(_mins,score[i]);
                _minp = std::min(_minp,pcore[i]);
            }
            eid_ming[e_id] = _mins;
            eid_minh[e_id] = _minp;
        
            if (log){
                std::cout<< "Edge ";
                for(auto u: elem.second)  std::cout<< u<<" ";
                std::cout<<": "<<_minp<<","<<_mins<<"\n";
            }
        }
        intvec temp(N); 
        for(size_t i = 0; i<N; i++){
            intvec vals;
            for (auto e_id : inc_edges[i]){
                if(eid_minh[e_id]>= pcore[i]) vals.push_back(eid_ming[e_id]);
            }
            temp[i] = hIndex(vals);
            if (temp[i]!=score[i]) flag = false;
        }
        for(size_t i = 0; i<N; i++) {
            score[i] = temp[i];
        }
        if (flag)   break;
     }
    a.output["total iteration(s)"] = std::to_string(iterations);

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        a.secondcore[node] = score[i];
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
}
