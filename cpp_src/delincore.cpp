// #include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <ctime>
#include <set>
#include <random>
#include <cassert>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"

typedef std::map<std::string, std::string> strstrMap;
template <typename T>
bool isSubsetOrEqual(std::set<T> const& a, std::set<T> const& b) {
   for(auto const& av:a){
      if(std::find(b.begin(),b.end(),av)==b.end())
          return false;
   }
   return true;
}
// int check_condition(Hypergraph& h, intIntMap& core){
//     /*
//     Checks that the sub-hypergraph induced by all nodes v with c(v)>=k has at least k incident hyperedges
//     in that sub-hypergraph \forall k \in [min_v c(v), max_v c(v) ]. (Coreness condition)
//     */
//     std::set<size_t> core_values;
//     int incorrect = 0;
//     for (auto node : h.init_nodes){
//         core_values.insert(core[node]);
//     }
//     bool condition_true = true;
//     for(auto p: core_values){
//         std::set<std::string> subnodes;
//         for(auto node: h.init_nodes){
//             auto node_str = std::to_string(node);
//             if (core[node]>= p){
//                 subnodes.insert(node_str);
//             }
//         }
//         Hypergraph subH;
//         size_t count = 0;
//         for(auto y: h.hyperedges){
//             std::vector<std::string> strvecE(y.size());
//             for(auto u:y){
//                 strvecE.push_back(std::to_string(u));
//             }
//             if (isSubsetOrEqual(std::set<std::string>(strvecE.begin(),strvecE.end()),subnodes)){
//                 subH.addEdge(count++,strvecE);
//             }
//         }
//         intIntMap node_deg; 
//         for(auto y: subH.hyperedges){
//             for(auto u:y){
//                 if (node_deg.find(u) == node_deg.end())
//                     node_deg[u] = 1;
//                 else 
//                     node_deg[u] += 1;
//             }
//         }
//         condition_true = true;
//         std::string violoating_node;
//         for(auto node: subH.init_nodes){
//             if (node_deg[node]< p){
//                 condition_true = false;
//                 violoating_node = node;
//             }
//         }
//         if (!condition_true){
//             std::cout<< "violating core: ("<<p<<")\n";
//             std::cout << "violating node: "<<violoating_node<<"\n";
//             incorrect+=1;
//         }
//     }
//     if (condition_true){
//         std::cout<< "No violation\n";
//     }
//     return incorrect;
// }

void residualhypergraph(intvec nodestodel, Hypergraph& h0, Hypergraph& h1){
    // For every edge e_id incident on v, remove e_id from every vertex u in edge[e_id] distinct from v
    std::set<std::string> subnodes; 
    for(auto u: nodestodel)    subnodes.insert(std::to_string(u));
    size_t count = 0;
    for(auto y: h0.hyperedges){
        std::vector<std::string> strvecE;
        bool add = true;
        for(auto u:y){
            if (subnodes.find(std::to_string(u)) == subnodes.end())
                strvecE.push_back(std::to_string(u));
            else{
                add = false;
                break;
            }
        }
        if (add)
            h1.addEdge(count++,strvecE);
    }
    h1.initialise();
}

void extract_nodes_to_delete(Hypergraph& h, intIntMap &core, int to_del, intvec& nodesto_del){
    /* If to_del == -1 deletes the entire innermost core,
    other-wise delete to_del number of nodes from innercore
    */
    srand(time(0));
    std::cout<<"deleting.\n";
    size_t max_core = std::numeric_limits<size_t>::min();
    for(auto pr: core){
        auto v = pr.first;
        size_t c = pr.second;
        max_core = std::max(c,max_core);
    }
    for(auto pr: core){
        if (pr.second == max_core)  nodesto_del.push_back(pr.first);
    }
    if (to_del > 0){
        assert (nodesto_del.size()>= to_del);
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(nodesto_del.begin(), nodesto_del.end(), g);
        if(nodesto_del.size()>=to_del){
            nodesto_del.resize(to_del);
        }
    }
    std::cout<<"#deleted: "<<nodesto_del.size()<<"\n";
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
            Hypergraph h1;
            std::string name = "";
            if (argc>=3){
                getHg(argv[2],h);
                name = argv[2];
                h.dataset = name+"_h0";
            }
            std::string init_type = "nbr"; // or "lub" (local upper bound)
            h.initialise();
            std::string alg;
            if (argc>=4){
                alg = argv[3];
            }
            int to_del;
            intvec nodesto_del;
            if(argc>=5){
                to_del = atoi(argv[4]);
            }
            else{
                to_del = -1; // -1 means delete the whole innermost core.
            }
            
            std::cout << argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
            if (argc>=6){
                // std::cout<<"Writing before-delete nbr: \n";
                if (atoi(argv[5])!=0)   h.writeneighborhood("../python_src/sirdata/log_"+h.dataset+"_"+std::to_string(to_del)+".csv");
            }
            if (alg == "Local-core-OPTIV" || alg == "naive_nbr"){    
                std::cout <<"Local-core-OPTIV (CSR + Modify core number at each step + bound + core correct) \n";
                Algorithm a(h);
                local_core_OPTIV(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                a.writecore("../python_src/sirdata/");
                extract_nodes_to_delete(h,a.core,to_del, nodesto_del);
                h1.dataset = name + "_h1";
                // std::cout<<"hypergraph before del\n";
                // h.printHypergraph();
                residualhypergraph(nodesto_del, h, h1);
                Algorithm a1(h1);
                local_core_OPTIV(h1.dataset, h1.hyperedges, h1.init_nodes, h1.node_index, a1, log);
                a1.writecore("../python_src/sirdata/");
                // std::cout<<"residual hypergraph: \n";
                // h1.printHypergraph();
            }
            if (alg == "deg"){    
                std::cout <<"naive_deg\n";
                Algorithm a(h);
                degreePeel(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                // a.output["num_threads"] = std::to_string(num_threads);
                // a.output["total iteration"] = std::to_string(0);
                // std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                a.writecore("../python_src/sirdata/");
                extract_nodes_to_delete(h,a.core,to_del, nodesto_del);
                h1.dataset = name + "_h1";
                // std::cout<<"hypergraph before del\n";
                // h.printHypergraph();
                residualhypergraph(nodesto_del, h, h1);
                Algorithm a1(h1);
                degreePeel(h1.dataset, h1.hyperedges, h1.init_nodes, h1.node_index, a1, log);
                a1.writecore("../python_src/sirdata/");
                // std::cout<<"residual hypergraph: \n";
                // h1.printHypergraph();
            }

            if(alg == "clique"){
                getClique(h);
                h.initialise();
                // h.printHypergraph();
                // auto ck_time = double(clock()-ck_start)/double(CLOCKS_PER_SEC);
                Algorithm a(h);
                local_core_clique(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, log);
                a.writecore("../python_src/sirdata/");
                extract_nodes_to_delete(h,a.core,to_del,nodesto_del);
                h1.dataset = name + "_h1";
                // std::cout<<"hypergraph before del\n";
                // h.printHypergraph();
                residualhypergraph(nodesto_del, h, h1);
                Algorithm a1(h1);
                local_core_clique(h1.dataset, h1.hyperedges, h1.init_nodes, h1.node_index, a1, log);
                a1.writecore("../python_src/sirdata/");
                // std::cout<<"residual hypergraph: \n";
                // h1.printHypergraph();
                // check_condition( h, a.core);
                // a.output["algo"] = alg;
                // a.output["num_threads"] = std::to_string(num_threads);
                // a.output["total iteration"] = std::to_string(1);
                // // a.output["execution time"] += ck_time;
                // a.output["init_time"] += ck_time;
                // a.write_results();
            }

            if (argc>=6){
                // std::cout<<"Writing before-delete nbr: \n";
                if (atoi(argv[5])!=0)   h1.writeneighborhood("../python_src/sirdata/log_"+h1.dataset+"_"+std::to_string(to_del)+".csv");
            }
        }
    }
}
