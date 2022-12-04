// #include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <ctime>
#include <set>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"

// typedef std::map<std::string, std::string> strstrMap;
typedef std::map<std::string, std::string> strstrMap;
template <typename T>
bool isSubsetOrEqual(std::set<T> const& a, std::set<T> const& b) {
   for(auto const& av:a){
      if(std::find(b.begin(),b.end(),av)==b.end())
          return false;
   }
   return true;
}
int check_conditiondeg(Hypergraph& h, intIntMap& core){
    /*
    Checks that the sub-hypergraph induced by all nodes v with c(v)>=k has at least k incident hyperedges
    in that sub-hypergraph \forall k \in [min_v c(v), max_v c(v) ]. (Coreness condition)
    */
    std::set<size_t> core_values;
    int incorrect = 0;
    for (auto node : h.init_nodes){
        core_values.insert(core[node]);
    }
    bool condition_true = true;
    for(auto p: core_values){
        std::set<std::string> subnodes;
        for(auto node: h.init_nodes){
            auto node_str = std::to_string(node);
            if (core[node]>= p){
                subnodes.insert(node_str);
            }
        }
        Hypergraph subH;
        size_t count = 0;
        for(auto y: h.hyperedges){
            std::vector<std::string> strvecE(y.size());
            for(auto u:y){
                strvecE.push_back(std::to_string(u));
            }
            if (isSubsetOrEqual(std::set<std::string>(strvecE.begin(),strvecE.end()),subnodes)){
                subH.addEdge(count++,strvecE);
            }
        }
        intIntMap node_deg; 
        for(auto y: subH.hyperedges){
            for(auto u:y){
                if (node_deg.find(u) == node_deg.end())
                    node_deg[u] = 1;
                else 
                    node_deg[u] += 1;
            }
        }
        condition_true = true;
        std::string violoating_node;
        for(auto node: subH.init_nodes){
            if (node_deg[node]< p){
                condition_true = false;
                violoating_node = node;
            }
        }
        if (!condition_true){
            std::cout<< "violating core: ("<<p<<")\n";
            std::cout << "violating node: "<<violoating_node<<"\n";
            incorrect+=1;
        }
    }
    if (condition_true){
        std::cout<< "No violation\n";
    }
    return incorrect;
}
int check_conditionnbr(Hypergraph& h, intIntMap& core){
    /*
    Checks that the sub-hypergraph induced by all nodes v with c(v)>=k has at least k neighhbors
    in that sub-hypergraph \forall k \in [min_v c(v), max_v c(v) ]. (Coreness condition)
    */
    std::set<size_t> core_values;
    int incorrect = 0;
    for (auto node : h.init_nodes){
        core_values.insert(core[node]);
    }
    bool condition_true = true;
    for(auto p: core_values){
        // auto p = x.second;
        // std::cout<<p<<" - "<<s<<"\n";
        std::set<std::string> subnodes;
        for(auto node: h.init_nodes){
            // if (core[node]>= p){
            if (core[node]>= p){
                auto node_str = std::to_string(node);
                subnodes.insert(node_str);
            }
        }
        // for(auto u: subnodes)   std::cout<<u<<" ";
        // std::cout<<"\n";
        Hypergraph subH;
        size_t count = 0;
        for(auto y: h.hyperedges){
            std::vector<std::string> strvecE(y.size());
            for(auto u:y){
                strvecE.push_back(std::to_string(u));
            }
            if (isSubsetOrEqual(std::set<std::string>(strvecE.begin(),strvecE.end()),subnodes)){
                subH.addEdge(count++,strvecE);
            }
        }
        // subH.initialise();
        // count neighbors
        intuSetintMap init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
        // intintvec edges( e_id_to_edge.size() ,intvec{});
        for(auto elem:subH.hyperedges){
            // auto elem = e_id_to_edge[eid];
            auto edge_sz = elem.size();
            for(auto v_id: elem){
                if(init_nbr.find(v_id) == init_nbr.end()){
                    init_nbr[v_id] = std::unordered_set<size_t>();
                }
                else{
                    auto _tmp = &init_nbr[v_id];
                    for (auto u: elem){
                        if (u!=v_id){
                            _tmp->insert(u);
                        }
                    }
                }
            }
        }
        condition_true = true;
        std::string violoating_node;
        for(auto node: subH.init_nodes){
            // if (subH.init_nbrsize[node]< p){
            if (init_nbr[node].size() < p){
                condition_true = false;
                violoating_node = node;
            }
        }
        if (!condition_true){
            std::cout<< "violating core: ("<<p<<")\n";
            std::cout << "violating node: "<<violoating_node<<"\n";
            incorrect+=1;
            // std::cout<<"sub-hyp.\n";
            // subH.printHypergraph();
            // std::cout<<"violating node nbr: ";
            // for (auto u: subH.init_nbr[violoating_node]) std::cout<< u<<" ";
            // std::cout<<"\n";
        }
    }
    if (condition_true){
        std::cout<< "No violation\n";
    }
    return incorrect;
}
int check_conditionkd(Hypergraph& h, intIntMap& core, intIntMap& secondcore){
    /*
    Checks that the sub-hypergraph induced by all nodes v with c(v)>=k has at least k neighhbors
    in that sub-hypergraph \forall k \in [min_v c(v), max_v c(v) ]. (Coreness condition)
    */
    std::set<std::pair<size_t,size_t>> core_pairs;
    int incorrect = 0;
    for (auto node : h.init_nodes){
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
                auto node_str = std::to_string(node);
                subnodes.insert(node_str);
            }
        }
        Hypergraph subH;
        size_t count = 0;
        for(auto y: h.hyperedges){
            std::vector<std::string> strvecE(y.size());
            for(auto u:y){
                strvecE.push_back(std::to_string(u));
            }
            if (isSubsetOrEqual(std::set<std::string>(strvecE.begin(),strvecE.end()),subnodes)){
                subH.addEdge(count++,strvecE);
            }
        }
        intuSetintMap init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
        std::map<size_t,size_t> node_deg;
        for(auto elem:subH.hyperedges){
            // auto elem = e_id_to_edge[eid];
            auto edge_sz = elem.size();
            for(auto v_id: elem){
                if(init_nbr.find(v_id) == init_nbr.end()){
                    init_nbr[v_id] = std::unordered_set<size_t>();
                }
                else{
                    auto _tmp = &init_nbr[v_id];
                    for (auto u: elem){
                        if (u!=v_id){
                            _tmp->insert(u);
                        }
                    }
                }
                if (node_deg.find(v_id) == node_deg.end())
                    node_deg[v_id] = 1;
                else 
                    node_deg[v_id] += 1;
            }
        }
        size_t violoating_node;
        for(auto node: subH.init_nodes){
            // if (subH.init_nbrsize[node]< p){
            if (init_nbr[node].size() < p || node_deg[node]< s){
                condition_true = false;
                violoating_node = node;
            }
        }
        if (!condition_true){
            std::cout<< "violating core: ("<<p<<","<<s<<")\n";
            std::cout << "violating node: "<<violoating_node<<"\n";
            incorrect+=1;
            // std::cout<<"sub-hyp.\n";
            // subH.printHypergraph();
            // std::cout<<"violating node nbr: ";
            // for (auto u: subH.init_nbr[violoating_node]) std::cout<< u<<" ";
            // std::cout<<"\n";
        }
    }

    // for(auto p: core_values){
    //     // auto p = x.second;
    //     // std::cout<<p<<" - "<<s<<"\n";
    //     std::set<std::string> subnodes;
    //     for(auto node: h.init_nodes){
    //         // if (core[node]>= p){
    //         if (core[node]>= p){
    //             auto node_str = std::to_string(node);
    //             subnodes.insert(node_str);
    //         }
    //     }
    //     // for(auto u: subnodes)   std::cout<<u<<" ";
    //     // std::cout<<"\n";
    //     Hypergraph subH;
    //     size_t count = 0;
    //     for(auto y: h.hyperedges){
    //         std::vector<std::string> strvecE(y.size());
    //         for(auto u:y){
    //             strvecE.push_back(std::to_string(u));
    //         }
    //         if (isSubsetOrEqual(std::set<std::string>(strvecE.begin(),strvecE.end()),subnodes)){
    //             subH.addEdge(count++,strvecE);
    //         }
    //     }
    //     // subH.initialise();
    //     // count neighbors
    //     intuSetintMap init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    //     // intintvec edges( e_id_to_edge.size() ,intvec{});
    //     for(auto elem:subH.hyperedges){
    //         // auto elem = e_id_to_edge[eid];
    //         auto edge_sz = elem.size();
    //         for(auto v_id: elem){
    //             if(init_nbr.find(v_id) == init_nbr.end()){
    //                 init_nbr[v_id] = std::unordered_set<size_t>();
    //             }
    //             else{
    //                 auto _tmp = &init_nbr[v_id];
    //                 for (auto u: elem){
    //                     if (u!=v_id){
    //                         _tmp->insert(u);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     condition_true = true;
    //     std::string violoating_node;
    //     for(auto node: subH.init_nodes){
    //         // if (subH.init_nbrsize[node]< p){
    //         if (init_nbr[node].size() < p){
    //             condition_true = false;
    //             violoating_node = node;
    //         }
    //     }
    //     if (!condition_true){
    //         std::cout<< "violating core: ("<<p<<")\n";
    //         std::cout << "violating node: "<<violoating_node<<"\n";
    //         incorrect+=1;
    //         // std::cout<<"sub-hyp.\n";
    //         // subH.printHypergraph();
    //         // std::cout<<"violating node nbr: ";
    //         // for (auto u: subH.init_nbr[violoating_node]) std::cout<< u<<" ";
    //         // std::cout<<"\n";
    //     }
    // }
    if (condition_true){
        std::cout<< "No violation\n";
    }
    return incorrect;
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
            if (argc>=3){
                getHg(argv[2],h);
                h.dataset = argv[2];
            }
            std::string init_type = "nbr"; // or "lub" (local upper bound)
            h.initialise();
            std::string alg;
            int iterations;
            bool log = false;
            if (argc>=4){
                alg = argv[3];
            }
            else{
               alg = "Peel";
            }
            if(argc>=5){
                iterations = atoi(argv[4]);
            }
            else
            iterations = 1;
            if(argc>=6){
                std::string s = argv[5];
                if(s[0]=='1') 
                log = true;
            }
            std::cout << argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
            for(int i=1;i<=iterations;i++){
                if (alg == "Local-core"){    //This is the optimised algorithm with opt_local_core_correct
                    std::cout <<"Sequential Local Core \n";
                    Algorithm a(h);
                    local_core(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(num_threads);
                    // a.printcore();
                    a.write_results();
                    // a.writecore();
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    if(i==iterations && log){
                        a.writecore();
                        a.writelog();
                    }
                }
                if (alg == "Local-core-OPTI"){    //This is the optimised algorithm with opt_local_core_correct
                    std::cout <<"Local-core-OPTI (CSR) \n";
                    Algorithm a(h);
                    local_core_OPTI(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(num_threads);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    a.write_results();
                   if(i==iterations && log){
                        a.writecore();
                        a.writelog();
                    }
                }
                if (alg == "Local-core-OPTII"){    //This is the optimised algorithm with opt_local_core_correct
                    std::cout <<"Local-core-OPTII (CSR + Modify core number at each step) \n";
                    Algorithm a(h);
                    local_core_OPTII(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(num_threads);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    a.write_results();
                    if(i==iterations && log){
                        a.writecore();
                        a.writelog();
                    }
                }
                if (alg == "Local-core-OPTIII"){    //This is the optimised algorithm with opt_local_core_correct
                    std::cout <<"Local-core-OPTIII (CSR + Modify core number at each step + bound) \n";
                    Algorithm a(h);
                    local_core_OPTIII(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(num_threads);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    a.write_results();
                    if(i==iterations && log){
                        a.writecore();
                        a.writelog();
                    }
                }
                if (alg == "Local-core-OPTIV"){    
                    std::cout <<"Local-core-OPTIV (CSR + Modify core number at each step + bound + core correct) \n";
                    Algorithm a(h);
                    local_core_OPTIV(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(num_threads);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    a.write_results();
                    if(i==iterations && log){
                        std::cout<<"Sum(Nu - cu) "<<a.nu_cu<<"\n";
                        a.writecore();
                        a.writelog();
                    }
                }
                if (alg == "Peel"){    
                    std::cout <<"Peel\n";
                    Algorithm a(h);
                    Peel(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(num_threads);
                    a.output["total iteration"] = std::to_string(0);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    a.write_results();
                    if(i==iterations && log){
                        std::cout<<"Execution time "<< a.exec_time<<"\n";
                        a.writecore();
                        a.writeNbrQ();
                    }
                }
                if (alg == "E-Peel"){    
                    std::cout <<"E-Peel\n";
                    Algorithm a(h);
                    EPeel(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(num_threads);
                    a.output["total iteration"] = std::to_string(0);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    a.write_results();
                    // a.writecore();
                    if(i==iterations && log){
                        std::cout<<"Execution time "<< a.exec_time<<"\n";
                        a.writecore();
                        a.writeNbrQ();
                    }
                }
                if (alg == "deg"){    
                    // if (i>1) continue;
                    std::cout <<"naive_degree\n";
                    Algorithm a(h);
                    degreePeel(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(num_threads);
                    a.output["total iteration"] = std::to_string(0);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    a.write_results();
                    if(i==iterations && log){
                        a.writecore();
                        check_conditiondeg( h, a.core);
                    }
                }
                if (alg == "clique"){
                    // if (i>1) continue;
                    std::cout <<"clique_core\n";
                    clock_t ck_start = clock();
                    getClique(h);
                    h.initialise();
                    auto ck_time = double(clock()-ck_start)/double(CLOCKS_PER_SEC);
                    Algorithm a(h);
                    local_core_clique(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.output["num_threads"] = std::to_string(1);
                    double totalinit_tm = atof(a.output["init_time"].c_str()) + ck_time;
                    a.output["init_time"] = std::to_string(totalinit_tm);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                    a.write_results();
                    if(i==iterations && log){
                        a.writecore();
                        check_conditionnbr( h, a.core);
                    }
                }
                if (alg=="kdcore"){
                    std::cout<< "Hybrid variant of k,d-core\n";
                    Algorithm a(h);
                    kdCorehybrid(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                    a.writekdcore();
                    check_conditionkd(h,a.core,a.secondcore);
                }
            }
        }
    }
    
}
