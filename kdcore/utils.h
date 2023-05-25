#ifndef UTILS_H
#define UTILS_H
#include <cstddef>
#include <cassert>
#include "hypergraph.h"

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
size_t hIndex_csr(size_t l, size_t r, size_t nbrs_f[], intvec & pcore) {
        if(l==r)
            return 0;
        size_t n = r-l;
        intvec hash(n + 1, 0);
        for(size_t i = l; i < r; ++i){
            size_t k = pcore[nbrs_f[i]];
            if(k >= n)
                hash[n]++;  
            else
                hash[k]++;
        }
        size_t paper = 0;
        for(size_t i = n; i >= 0; --i){
            paper += hash[i];
            if(paper >= i)
                return i;
        }
        return -1;
    }
#endif