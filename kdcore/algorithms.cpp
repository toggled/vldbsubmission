#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include <algorithm>
#include <cassert>
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
void Algorithm::writecore(std::string folder="../output/"){
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
    for(auto i=output.begin(); i != output.end();++i )
        std::cout<< i->first<<",";
    
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
void local_kdcore( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a){
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
    
    start = clock();
    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrsizes[i]; // initialize pcore
        score[i] =  inc_edges[i].size();
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
        if (flag)
            break;
    
    }
    a.output["total iteration(p)"] = std::to_string(iterations);

    // Compute secondary core-numbers
    iterations = 0;
    intIntMap eid_minh;
    intIntMap eid_ming;

    // intIntMap monotonicity;
     while (1){
        iterations+=1;
        // std::cout<<iterations<<"\n";
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
        }
        intvec temp(N); 
        for(size_t i = 0; i<N; i++){
            // size_t degSubg=0;
            intvec vals;
            for (auto e_id : inc_edges[i]){
                // if(eid_minh[e_id]>= pcore[i])   degSubg +=1;
                if(eid_minh[e_id]>= pcore[i]) vals.push_back(eid_ming[e_id]);
            }
            // score[i] = degSubg;
            auto Hvalue = hIndex(vals);
            
            // if (Hvalue<= score[i]) {    
            //         if (iterations == 1) monotonicity[i] = true;
            //         else{
            //             if (!monotonicity[i]) {}
            //         }
            // }else{
            //     monotonicity[i] = false;
            // }
            // // size_t count = 0;
            // // for (auto e_id : inc_edges[i]){ 
            // //     if (eid_ming[e_id]>= Hvalue)    count++;
            // // }
            // // if (count != score[i]){
            // if (Hvalue != score[i]){
            //     score[i] = Hvalue,; 
            //     // score[i] = count; 
            //     flag = false;
            // }
            // temp[i] = std::min(Hvalue,score[i]);
            temp[i] = Hvalue;
            if (temp[i]!=score[i]) flag = false;
        }
        for(size_t i = 0; i<N; i++) {
            score[i] = temp[i];
        }
        // if (log){
        //      for(size_t i = 0; i<N; i++)
        //         std::cout<< init_nodes[i]<<":"<< score[i]<<"\n";
        // }

        // for(size_t i = 0; i<N; i++){
        //     intvec vals(inc_edges[i].size());
        //     size_t ii=0;
        //     for (auto e_id : inc_edges[i]){
        //         vals[ii++] = eid_minh[e_id];
        //     }
        //     size_t H_value = hIndex(vals);
        //     if (H_value < score[i]){
        //         score[i] = H_value;
        //         flag = false;
        //     }
        //     // else gn[i] = score[i];
        //     std::cout<< init_nodes[i]<<":"<< score[i]<<"\n";
        // }
        if (flag)
            break;
     }
    a.output["total iteration(s)"] = std::to_string(iterations);
    end_main = clock();

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        a.secondcore[node] = score[i];
        // std::cout << i<< ": "<< node << "->"<< pcore[i] <<" , "<<score[i]<<"\n";
    }
    // std::cout<<"monotonicity: \n";
    // for(size_t i = 0; i<N; i++) std::cout<<monotonicity[i]<<" ";
    // std::cout<<"\n";
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
}
void local_kdcoreTest( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a){
    std::cout<<"start of local_kdcore\n";
    a.output["algo"] = "kdcore";
    clock_t start, end;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t N = init_nodes.size();
    intvec pcore(N); //
    intvec score(N); //
    std::set<size_t> A; // Initialization for gaberts algorithm.
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    strInthashMap node_index; // (use hashtable instead of dictionary => Faster on large |V| datasets.)
    for(size_t i = 0; i<N; i++) {
        node_index[init_nodes[i]] = i; // initialize node_index
        A.insert(i); /* Initializatoin for gaberts Alg.*/
    }
    intvec nbrsizes(N); //
    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    std::vector<intvec> nbrs(N);
    // compute initial neighbors and number of neighbors
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    // compute initial neighbors and number of neighbors
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

    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    nbrs_N[0] = 0;
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));

    size_t glb = std::numeric_limits<size_t>::max();
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
    }


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
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();

    // initialise core to a upper bound
    // #pragma omp parallel for default(none) shared(nbrsizes,pcore,N,glb,llb)
    for (size_t i = 0; i < N; i++){
        // printf("%d processing: %d\n",omp_get_thread_num(),i);
        pcore[i] = nbrsizes[i]; // initialize pcore
        llb[i] = std::max(llb[i],glb);
        score[i] = inc_edges_N[i+1]-inc_edges_N[i];
        assert (score[i]==inc_edges[i].size());
    }
        // intvec hn(N);
    intvec gn(N);
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
        // for every edge update is minimum of hindex(constituent vertices)
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

    // while (1){
    //     iterations+=1;
    //     // std::cout<<"iteration: "<<iteration<<"\n";
    //     // for(size_t i=0; i<N; i++)
    //     // {
    //     //     auto node = hg.init_nodes[i];
    //     //     std::cout << i<< ": "<< node << "->"<< pcore[i] <<"\n";
    //     // }
        
    //     bool flag = true;
    //     // compute h-index and update core
    //     for(size_t i = 0; i<N; i++){
    //         intvec vals(nbrsizes[i]);
    //         size_t ii=0;
    //         for (auto u_id : nbrs[i]){
    //             vals[ii++] = pcore[u_id];
    //         }
    //         size_t H_value = hIndex(vals);
    //         if (H_value < pcore[i]) 
    //         hn[i] = H_value;
    //         else hn[i] = pcore[i];
    //     }
        
    //     for (size_t i = 0; i<N; i++){
    //         bool lccsat = LCCSAT_check(inc_edges[i],edges,i,hn[i],hn);
    //         if (lccsat == false){ 
    //             flag = false;   //Why this step
    //             auto hhatn = core_correct(inc_edges[i],edges,i,hn[i],hn);
    //             pcore[i]  = hhatn;
    //         }else{
    //             pcore[i] = hn[i];
    //         }
    //     }
    //     // time_t end1 = clock();
    //     // if (log){
    //     //     strstrMap h0;
    //     //     for(int i=0;i<N;i++){
    //     //         h0[init_nodes[i]] = std::to_string(pcore[i]);
    //     //     }
    //     //     h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
    //     //     a.hnlog.push_back(h0);
    //     // }
    //     if (flag)
    //         break;
    
    // }

    a.output["total iteration(p)"] = std::to_string(iterations);
    // for(size_t i = 0; i< M; i++){
    //     size_t _min = N+1;      //Why M+1 and why calculate 
    //     for (auto u_id: edges[i])   _min = std::min(_min,pcore[u_id]);
    //     min_e_hindex[i] = _min;
    // }
    // Compute secondary core-numbers
    iterations = 0;
    intIntMap eid_minh;
    intIntMap eid_ming;
    // if (log){
    //     std::cout<<"primary core: \n";
    //     for(size_t i=0; i<N; i++)
    //     {
    //         std::cout<<init_nodes[i]<<": "<<pcore[i]<<"\n";
    //     }
    // }
    bool log = true;
    // Secondary-core computation loop
    while (1){
        iterations+=1;
        bool flag = true;
        for (auto elem: e_id_to_edge){
            auto e_id = elem.first;
            auto _tmp = intvec();
            size_t _mins = std::numeric_limits<size_t>::max();
            size_t _minp = std::numeric_limits<size_t>::max();
            // bool same_pcore = true;
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
            if (temp[i]!=score[i]) {
                if (log){   std::cout<< init_nodes[i]<<" "<<score[i]<<" => "<<temp[i]<<"\n"; }
                flag = false;
            }
        }
        for(size_t i = 0; i<N; i++) {
            score[i] = temp[i];
        }
        if (flag)   break;
     }

    // /* Apply gaberts algorithm. */
    // while (A.size()){
    //     iterations+=1;
    //     // compute h-index and update core
    //     std::set<size_t> Ap;
    //     for(auto i: A){
    //         bool local_flag = false;
    //         intvec vals;
    //         for(auto e_id: inc_edges[i]){
    //             if (min_e_hindex[e_id]>=pcore[i]){
    //                 size_t min_h = std::numeric_limits<size_t>::max();
    //                 for(auto u: edges[e_id]){
    //                     if (u!=i){
    //                         min_h = std::min(min_h,score[u]);
    //                     }    
    //                 }
    //                 vals.push_back(min_h);
    //             }
    //         }
    //         auto tmp = hIndex(vals);
    //         if (tmp!=score[i])  {
    //             local_flag = true;
    //         }
    //         score[i] = tmp;
    //         if (local_flag){
    //             Ap.insert(i);
    //             // Insert all neighors j of i such that
    //             // pcore[j] >= pcore[i] and edge containing i,j has min-p hindex >= pcore[i]
    //             for(auto e_id: inc_edges[i]){
    //                 if (min_e_hindex[e_id]>=pcore[i]){
    //                     for(auto u:  edges[e_id]){
    //                         if (u!=i)   Ap.insert(u);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     A = std::move(Ap);
    // }

    a.output["total iteration(s)"] = std::to_string(iterations);

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        a.secondcore[node] = score[i];
        std::cout <<node<<","<<a.core[node]<<","<<a.secondcore[node]<<".\t";
    }
    std::cout<<"\n";
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
}
