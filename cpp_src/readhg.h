#ifndef READHG_H
#define READHG_H
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <iostream>
#include <iterator>
#include <vector>
#include "hypergraph.h"

std::map <std::string,std::string> dataset_to_filename = {
            {"default" , "../data/datasets/synthetic/default.hyp"},
            {"enron" , "../data/datasets/real/Enron.hyp"},
            {"congress" , "../data/datasets/real/congress-bills.hyp"},
            {"contact" , "../data/datasets/real/contact-primary-school.hyp"},
            {"dblp", "../data/datasets/real/DBLP.hyp"},
            {"protein", "../data/datasets/protein/humancomplexes.hyp"},
            {"aminer","../data/datasets/real/aminer.hyp"},
            {"klay","../data/datasets/kenneth_lay/klay_s.hyp"},           
           
            {"bin_2" , "../data/datasets/synthetic/binomial_5_500_4_0.200000_sample_2_iter_1.txt"},
            {"bin_5" , "../data/datasets/synthetic/binomial_5_500_3_0.200000_sample_5_iter_1.txt"},
            {"pref", "../data/datasets/synthetic/pref_1000000_3_1.hyp"}
        };

template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

void getHg(std::string dataname, Hypergraph & hg){
    std::cout<< dataname<<"\n";
    std::cout << dataset_to_filename[dataname]<<"\n";
    std::ifstream infile(dataset_to_filename[dataname]);
    std::string line;
    // std::map< int, std::vector<std::string>  > Edges;
    size_t i= 0;
    while (std::getline(infile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        // std::cout<< i<<" "<<line<<"\n";
        std::vector<std::string> x = split(line ,',');
        hg.addEdge(i,x);
        i++;
    }
    // hg.print();
}

void getClique(Hypergraph &hg){
    std::set<std::pair<std::string,std::string>> edges;
    for(auto e_id:hg.hyperedges){
        for(int i=0;i<e_id.size();i++){
            for(int j=i+1;j<e_id.size();j++){
                std::pair<size_t,size_t> edge = std::make_pair(std::min(e_id[i],e_id[j]),std::max(e_id[i],e_id[j]));
                std::pair<std::string,std::string> s_edge = std::make_pair(std::to_string(edge.first),std::to_string(edge.second));
                if(edges.find(s_edge)==edges.end()){
                    edges.insert(s_edge);
                }
            }
        }
    }
    hg.hyperedges.clear();
    hg.init_nodes.clear();
    hg.node_index.clear();
    size_t i = 0;
    for(auto edge : edges){
        std::vector<std::string> x = {edge.first,edge.second};
        // std::cout<<edge.first<<","<<edge.second<<"\n";
        hg.addEdge(i,x);
        i++;
    }
}
#endif