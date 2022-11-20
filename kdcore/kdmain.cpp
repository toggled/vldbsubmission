// #include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <ctime>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"

template <typename T>
bool isSubsetOrEqual(std::set<T> const& a, std::set<T> const& b) {
   for(auto const& av:a){
      if(std::find(b.begin(),b.end(),av)==b.end())
          return false;
   }
   return true;
}
int check_condition(Hypergraph& h, strIntMap& core, strIntMap& secondcore){
    std::set<std::pair<size_t,size_t>> core_pairs;
    int incorrect = 0;
    for (auto node : h.init_nodes){
        // std::cout<<node<<"\n";
        std::pair<size_t,size_t> x = std::make_pair(core[node],secondcore[node]);
        core_pairs.insert(x);
    }
    bool condition_true = true;
    for(auto x: core_pairs){
        auto p = x.first;
        auto s = x.second;
        std::set<std::string> subnodes;
        for(auto node: h.init_nodes){
            // if (core[node]>= p){
            if (core[node]>= p && secondcore[node]>= s){
                subnodes.insert(node);
            }
        }
        Hypergraph subH;
        size_t count = 0;
        for(auto y: h.e_id_to_edge){
            // auto id = y.first;
            auto strvecE = std::set<std::string>(y.second.begin(),y.second.end());
            if (isSubsetOrEqual(strvecE,subnodes)){
                subH.addEdge(count++,y.second);
            }
        }
        subH.initialise();
        subH.init_nbrs();
        condition_true = true;
        std::string violoating_node;
        for(auto node: subH.init_nodes){
            // if (subH.init_nbrsize[node]< p){
            if (subH.init_nbrsize[node]< p || subH.inc_dict[node].size() < s){
                condition_true = false;
                violoating_node = node;
            }
        }
        if (!condition_true){
            // std::cout<< "violating core: ("<<p<<","<<s<<")\n";
            // std::cout << "violating node: "<<violoating_node<<"\n";
            incorrect+=1;
        }
    }
    if (condition_true){
        std::cout<< "No violation\n";
    }
    else{
        // h.printHypergraph();
        // for (auto node : h.init_nodes){
        //     std::cout<<node<<": "<<"("<<core[node]<<","<<secondcore[node]<<")\n";
        // }
    }
    return incorrect;
}

int main(int argc, char *argv[])
{
    if (argc >= 2)
    {
            Hypergraph h;
            if (argc>=3){
                getHg(argv[2],h);
                h.dataset = argv[2];

            }
            std::string init_type = "nbr"; // or "lub" (local upper bound)
            h.initialise();
            std::string alg;
            int iterations=1;
            // bool log = false;
            if (argc>=4){
                alg = argv[3];
                std::cout<<argv[2]<<" - "<<argv[3]<<"\n";
            }
            else{
               alg = "Peel";
            }
            if(argc>=5){
                iterations = atoi(argv[4]);
            }
            // if(argc>=6){
            //     std::string s = argv[5];
            //     if(s[0]=='1') 
            //     log = true;
            // }
            std::cout <<"params: "<< argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
            std::cout<<"Iterations: "<<iterations<<"\n";
            for(int i=1;i<=iterations;i++){
                std::cout<<"i: "<<i<<"\n";
                if (alg == "kdcore"){  
                    std::cout<<"alg==kdcore\n";
                    Algorithm a(h);
                    local_kdcore(h.dataset, h.e_id_to_edge, h.inc_dict, h.init_nodes, a);
                    auto c = check_condition(h,a.core,a.secondcore);
                    std::cout<<"incorrect: "<<c<<"/"<<h.init_nodes.size()<<"\n";
                    // a.output["num_threads"] = std::to_string(num_threads);
                    // a.write_results();
                    // a.writecore();
                    // a.writelog();
                }
            }                
        }
}
    
