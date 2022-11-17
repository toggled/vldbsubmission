// #include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <ctime>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"

typedef std::map<std::string, std::string> strstrMap;

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
            int iterations;
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
            else
            iterations = 1;
            // if(argc>=6){
            //     std::string s = argv[5];
            //     if(s[0]=='1') 
            //     log = true;
            // }
            std::cout << argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
            std::cout<<"Iter: "<<iterations<<"\n";
            for(int i=1;i<=iterations;i++){
                if (alg == "kdcore"){  
                    std::cout<<"alg==kdcore\n";
                    Algorithm a(h);
                    local_kdcore(h.dataset, h.e_id_to_edge, h.inc_dict, h.init_nodes, a);
                    // a.output["num_threads"] = std::to_string(num_threads);
                    a.write_results();
                    a.writecore();
                    // a.writelog();
                }
            }                
        }
}
    
