// #include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <ctime>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"

// typedef std::map<std::string, std::string> strstrMap;

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
                // if (alg == "E-Peel"){    
                //     std::cout <<"E-Peel\n";
                //     Algorithm a(h);
                //     EPeel(h.dataset, h.e_id_to_edge, h.inc_dict, h.init_nodes, a,log);
                //     a.output["num_threads"] = std::to_string(num_threads);
                //     a.output["total iteration"] = std::to_string(0);
                //     a.write_results();
                //     if(i==iterations && log){
                //         std::cout<<"Execution time "<< a.exec_time<<"\n";
                //         a.writecore();
                //         a.writeNbrQ();
                //     }
                // }

            }
        }
    }
    
}
