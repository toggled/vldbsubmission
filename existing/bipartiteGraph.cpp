#include <iostream>
#include <sstream>
#include <ctime>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"
// #include "utils.h"
// #include <omp.h>

typedef std::map<std::string, std::string> strstrMap;
size_t hIndex(intvec & citations) {
        if(citations.empty())
            return 0;
        size_t n = citations.size();
        intvec hash(n + 1, 0);
        for(size_t i = 0; i < n; ++i){
            if(citations[i] >= n)
                hash[n]++;  
            else
                hash[citations[i]]++;
        }
        size_t paper = 0;
        for(size_t i = n; i >= 0; --i){
            paper += hash[i];
            if(paper >= i)
                return i;
        }
        return -1;
    }
void getBipartite(Hypergraph &hg, strvIntMap &v_to_eid, intsStrMap &eid_to_v){

    for(auto e_id:hg.e_id_to_edge){
        
        for(auto v:e_id.second){
            v_to_eid[v].push_back(e_id.first);
            eid_to_v[e_id.first].insert(v);
        }
    }

}

std::set<std::string> getDist2nbr(strvIntMap &v_to_eid, intsStrMap &eid_to_v, std::string u){

    std::set<std::string> s;
    for(auto eid:v_to_eid[u]){
        for(auto v:eid_to_v[eid]){
            if(v!=u)
            s.insert(v);
        }
    }
    return s;

}

std::set<int> getDist2nbredge(strvIntMap &v_to_eid, intsStrMap &eid_to_v, int u){

    std::set<int> s;
    for(auto v:eid_to_v[u]){
        for(auto e:v_to_eid[v]){
            if(e!=u)
            s.insert(e);
        }
    }
    return s;

}

void bipartiteDeleteu(strvIntMap &v_to_eid, intsStrMap &eid_to_v, std::string u){

    for(auto eid:v_to_eid[u]){
        eid_to_v[eid].erase(u);
    }
    v_to_eid.erase(u);

}

void bipartitedist2core(Hypergraph &h){

    time_t start, end;
    start = clock();
    std::map<std::string,int> core;
    strvIntMap v_to_eid;
    intsStrMap eid_to_v;
    getBipartite(h,v_to_eid,eid_to_v);
    std::cout <<"E-Peel\n";
    std::map<int,std::set<std::string>> bucket;
    std::map<std::string,int> node_to_num_nbrs;
    int num_nodes = 0;

    for(int i=0;i<h.init_nodes.size();i++){
        std::string u = h.init_nodes[i];
        int len_neighbours = getDist2nbr(v_to_eid,eid_to_v,u).size();
        node_to_num_nbrs[u] = len_neighbours;
        bucket[len_neighbours].insert(u);
        num_nodes++;
    }

    for(int k=1;k<=num_nodes;k++){

        while(!bucket[k].empty()){
            auto it = bucket[k].end();
            it--;
            std::string v = *it;
            bucket[k].erase(it);
            std::set<std::string> nbr_v = getDist2nbr(v_to_eid,eid_to_v,v);
            bipartiteDeleteu(v_to_eid,eid_to_v,v);
            core[v] = k;
            for(auto u:nbr_v){

                int len_nbrs_u = getDist2nbr(v_to_eid,eid_to_v,u).size();
                int max_value = std::max(len_nbrs_u,k);
                int prev_idx = node_to_num_nbrs[u];
                bucket[prev_idx].erase(u);
                bucket[max_value].insert(u);
                node_to_num_nbrs[u] = max_value;

            }
        }

    }
    end = clock();

    std::cout<<"Execution time: "<<double(end - start) / double(CLOCKS_PER_SEC)<<"\n";

    std::string file = "output/core_E-Peel_bipartite_"+h.dataset+".csv";
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

void bipartitedist2core_local(Hypergraph &h, Algorithm &a){
    a.output["algo"] = "bipartite";
    std::cout<<"local-bp\n";
    // v is string
    time_t start, end;
    std::map<std::string,int> core;
    int num_edge = h.e_id_to_edge.size();
    strInthashMap node_index;
    int index = num_edge;
    for(int i=0;i<h.init_nodes.size();i++){
        node_index[h.init_nodes[i]] = index;
        index++;
    }
    start = clock();
    strvIntMap v_to_eid;
    intsStrMap eid_to_v;
    getBipartite(h,v_to_eid,eid_to_v);
    // std::cout <<"Local-core\n";
    intvec pcore(eid_to_v.size()+v_to_eid.size());

    //Initialise pcore
    for(auto v:h.init_nodes){
        pcore[node_index[v]] = getDist2nbr(v_to_eid,eid_to_v,v).size() + v_to_eid[v].size();
        // std::cout<<v<<" "<<pcore[node_index[v]]<<"\n";
    }
    // std::cout<<"edges\n";
    for(auto e:h.e_id_to_edge){
        pcore[e.first] = getDist2nbredge(v_to_eid,eid_to_v,e.first).size() + eid_to_v[e.first].size();
        // std::cout<<e.first<<" "<<pcore[e.first]<<"\n";
    }
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();
    int iterations = 0;
    while(1){
        // std::cout<<iterations<<"\n";
        iterations++;
        // for(int i=0;i<pcore.size();i++)
        // std::cout<<pcore[i]<<" ";
        // std::cout<<"\n";
        assert(eid_to_v.size()==h.e_id_to_edge.size());
        assert(v_to_eid.size()==h.init_nodes.size());
        intvec H(eid_to_v.size()+v_to_eid.size());
        for(auto v:h.init_nodes){

            //Given H(n-1), calculate A(n-1). 
            // attainability is per start node and end node with distance <= 2 from start node. 

            std::unordered_map<int,int> attainability;
            for(auto e:v_to_eid[v]){
                attainability[e] = pcore[e];
                for(auto u:eid_to_v[e]){
                    if(u!=v){
                        if(attainability.find(node_index[u])==attainability.end()){
                            attainability[node_index[u]] = std::min(pcore[e],pcore[node_index[u]]);
                        }else{
                            int curr_attain = std::min(pcore[e],pcore[node_index[u]]);
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

        for(auto e:h.e_id_to_edge){

            //Given H(n-1), calculate A(n-1). 
            // attainability is per start node and end node with distance <= 2 from start node. 

            std::unordered_map<int,int> attainability;
            for(auto v:eid_to_v[e.first]){
                attainability[node_index[v]] = pcore[node_index[v]];
                for(auto eid : v_to_eid[v]){
                    if(e.first!=eid){
                        if(attainability.find(eid)==attainability.end()){
                            attainability[eid] = std::min(pcore[eid],pcore[node_index[v]]);
                        }else{
                            int curr_attain = std::min(pcore[eid],pcore[node_index[v]]);
                            attainability[eid] = std::max(attainability[eid],curr_attain);
                        }
                    }
                }
            }

            // Use A(n-1) to calculate H(n).
            intvec vals;
            for(auto val : attainability)
            vals.push_back(val.second);
            H[e.first] = hIndex(vals);

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
    end = clock();
    a.output["total iteration"] = std::to_string(iterations);
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    std::cout<<"Execution time: "<<a.exec_time<<"\n";

    std::string file = "../output/core_bipartite_"+h.dataset+".csv";
    std::cout<<"writing core to file: "<<file<<"\n";
    std::stringstream ss;
    for(auto elem: h.init_nodes)
    {
        ss << elem << "," << pcore[node_index[elem]] << "\n";
    }
    std::ofstream out(file.c_str());
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();

    file = "../output/core_bipartite_"+h.dataset+"_edge.csv";
    std::stringstream ss2;
    for(auto elem: h.e_id_to_edge)
    {
        ss2 << elem.first << "," << pcore[elem.first] << "\n";
    }
    std::ofstream out2(file.c_str());
    if(out2.fail())
    {
        out2.close();
    }
    out2 << ss2.str();
    out2.close();

}

int main(int argc, char *argv[])
{
    if (argc >= 2)
    {
        std::istringstream iss( argv[1] );
        int num_threads;

        if (iss >> num_threads)
        {
            std::cout << num_threads<<"\n";
            Hypergraph h;
            // std::cout<<"argv[2]"<<argv[2]<<"\n";
            if (argc>=3){
                if(std::string(argv[2])=="rand"){
                    std::cout<<"Generating random hyp.\n";
                    // getrandomHg(h);
                    getrandomHg(h,8,5,5);
                    h.dataset = "rand";
                }
                else{
                    getHg(argv[2],h);
                    h.dataset = argv[2];
                }
            }
            else{
                getHg("default",h);
                h.dataset = "default";
            }
            std::string algo = "E-Peel";
            if(argc>=4){
                algo = argv[3];
            }
            std::cout<<"Bipartite graph\n";
            // std::string init_type = "nbr"; // or "lub" (local upper bound)
            h.initialise();
            if (argc>=5){
                if (atoi(argv[4])!=0)   h.writeneighborhood("../output/log_"+h.dataset+"_Nv.csv");
                return 0;
            }
            Algorithm a(h);
            // if(algo=="E-Peel"){
            //     bipartitedist2core(h);
            // }
            // else if(algo=="bipartite"){
            bipartitedist2core_local(h,a);
            // }
            a.output["num_threads"] = std::to_string(num_threads);
            // a.write_results();
            a.printcore();
        }
    }
}
// g++ -std=c++11 -g -o bpmain bipartiteGraph.cpp hypergraph.cpp  algorithms.cpp readhg.h