#ifndef READHG_H
#define READHG_H

#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <iostream>
#include <iterator>
#include <vector>
#include <cassert>
#include <random>
#include "hypergraph.h"

std::map <std::string,std::string> dataset_to_filename = {
            
            {"enron" , "../data/datasets/real/Enron.hyp"},
            {"congress" , "../data/datasets/real/congress-bills.hyp"},
            {"contact" , "../data/datasets/real/contact-primary-school.hyp"},
            {"dblp", "../data/datasets/real/DBLP.hyp"},
            {"protein", "../data/datasets/real/humancomplexes.hyp"},
            {"aminer","../data/datasets/real/aminer.hyp"},
            {"klay","../data/datasets/real/klay_s.hyp"},           
           
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
// See https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
std::vector<int> FisherYatesShuffle(std::size_t size, std::size_t max_size, std::mt19937& gen)
{
    assert(size <= max_size);
    std::vector<int> b(size);
 
    for(std::size_t i = 0; i != max_size; ++i) {
        std::uniform_int_distribution<> dis(0, i);
        std::size_t j = dis(gen);
        if(j < b.size()) {
            if (i < b.size()) {
                b[i] = b[j];
            }
            b[j] = i;
        }
    }
    return b;
}

void getrandomHg(Hypergraph & hg, int N = 7, int M = 10, int edge_size_ub=3, unsigned seed=1, bool log = false){
    size_t _numE = 0;
    std::vector<std::string> edges;
    if (M == 1){
        if (N <= edge_size_ub){
            // e = tuple([str(i) for i in range(1,n+1)])
            std::vector<std::string> e;
            for (int i=0; i<N; i++)    e.push_back(std::to_string(i));
            std::sort(e.begin(),e.end());
            std::string _tmp = "";
            for (auto i: e) _tmp+=(i+",");
            if (std::find(edges.begin(),edges.end(),_tmp) == edges.end()){
                edges.push_back(_tmp);
                if (log){
                    for(auto i:e)   std::cout<<i<<",";
                }
                hg.addEdge(_numE++, e);
                return;
            }
        }
    }
    // unsigned seed = unsigned(std::time(nullptr));
    std::srand(seed);
    std::random_device rd;
    // std::mt19937 gen(rd());
    std::mt19937 gen(seed);
    while (M){
        // std::cout <<M<<"-";
        std::uniform_int_distribution<int> uni(2,edge_size_ub); // Guaranteed unbiased
        size_t edge_sz = uni(gen);
        // std::cout<<"edge_sz: "<<edge_sz<<"\n";
        std::vector<int> e = FisherYatesShuffle(edge_sz, N, gen);
        std::vector<std::string> E;
        for (auto i: e)    E.push_back(std::to_string(i));
        std::sort(E.begin(),E.end());
        std::string _tmp = "";
        for (auto e: E) _tmp+=(e+",");
        // std::cout<<_tmp<<"\n";
        if (std::find(edges.begin(),edges.end(),_tmp) == edges.end()){
            edges.push_back(_tmp);
            if (log){
                    for(auto i:E)   std::cout<<i<<",";
                    std::cout<<"\n";
                }
            hg.addEdge(_numE++, E);
            M-=1;
        }
    }
}

#endif