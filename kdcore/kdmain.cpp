// #include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <ctime>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"

std::pair<size_t,size_t> GenhIndex(intvec & citations1,intvec & citations2) {
    if(citations1.empty() || citations1.empty())
        return std::make_pair(0,0);
    size_t n1 = citations1.size();
    size_t n2 = citations2.size();
    assert (n1==n2);
    intvec hash1(n1 + 1, 0);
    intvec hash2(n2 + 1, 0);
    for(size_t i = 0; i < n1; ++i){
        if(citations1[i] >= n1)
            hash1[n1]++;  
        else
            hash1[citations1[i]]++;
        
        if(citations2[i] >= n2)
            hash2[n2]++;  
        else
            hash2[citations2[i]]++;
    }
    size_t x = 0;
    size_t y = 0;
    size_t tmp1 = 0, tmp2 = 0;
    bool flag1 = false, flag2 = false;
    for(size_t i = n1; i >= 0; --i){
        if (flag1 && flag2) break;
        x += hash1[i];
        y += hash2[i];
        if(x >= i && !flag1){
            flag1 = true;
            tmp1 = i;
            tmp2 = y;  
            break;
        }
        if (y >= i && !flag2){
            flag2 = true;
            tmp2 = i;
            tmp1 = x;
            break;
        }
    }
    return std::make_pair(tmp1,tmp2);
}

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
        std::cout<<p<<" - "<<s<<"\n";
        std::set<std::string> subnodes;
        for(auto node: h.init_nodes){
            // if (core[node]>= p){
            if (core[node]>= p && secondcore[node]>= s){
                subnodes.insert(node);
            }
        }
        for(auto u: subnodes)   std::cout<<u<<" ";
        std::cout<<"\n";
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
            std::cout<< "violating core: ("<<p<<","<<s<<")\n";
            std::cout << "violating node: "<<violoating_node<<"\n";
            incorrect+=1;
            std::cout<<"sub-hyp.\n";
            subH.printHypergraph();
            std::cout<<"violating node nbr: ";
            for (auto u: subH.init_nbr[violoating_node]) std::cout<< u<<" ";
            std::cout<<"\n";
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
            intvec a ={2,2,3};
            intvec b = {1,1, 1};
            std::pair<size_t,size_t> out1 = GenhIndex(a,b);
            // std::cout <<"Gen-hindex: "<<out1.first<<" "<<out1.second<<"\n";
            // return 0;
            std::cout <<"params: "<< argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
            std::cout<<"Iterations: "<<iterations<<"\n";
            for(int i=1;i<=iterations;i++){
                if(std::string(argv[2])=="rand"){
                    std::cout<<"Generating random hyp.\n";
                    // getrandomHg(h);
                    getrandomHg(h,15,15,5,12);
                    h.dataset = "rand";
                    h.initialise();
                    h.printHypergraph();
                }

                std::cout<<"i: "<<i<<"\n";
                if (alg == "kdcore"){  
                    std::cout<<"alg==kdcore\n";
                    Algorithm a(h);
                    // local_kdcore(h.dataset, h.e_id_to_edge, h.inc_dict, h.init_nodes, a);
                    local_kdcoreTest(h.dataset, h.e_id_to_edge, h.inc_dict, h.init_nodes, a);
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
    
