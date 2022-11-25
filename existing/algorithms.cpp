#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include <algorithm>
#include "hypergraph.h"
#include "algorithms.h"
// #include "utils.h"

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
    std::string file = "output/core_"+output["algo"]+"_"+hg.dataset+".csv";
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
bool Algorithm::write_results(){
    std::string file = "../output/exresults.csv";
    std::stringstream ss;
    
    // // print the order in which the output keys are saved in the csv file.
    for(auto i=output.begin(); i != output.end();++i )
        std::cout<< i->first<<",";
    std::cout<<"\n";
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
    std::string file = "output/nbrq_results.csv";
    // std::cout<<file<<"\n";
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
    std::string file = "output/log_"+output["algo"]+"_"+hg.dataset+".csv";
    // std::cout<<file<<"\n";
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
