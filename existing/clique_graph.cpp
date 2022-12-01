#include <iostream>
#include <sstream>
#include <ctime>
#include "readhg.h"
#include "hypergraph.h"
#include "algorithms.h"
// #include "utils.h"

typedef std::map<std::string, std::string> strstrMap;

size_t hIndex_CSR(size_t l, size_t r, size_t nbrs_f[], intvec & pcore) {
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

void local_core_clique( std::string dataset, std::map<size_t, strvec > &e_id_to_edge, std::map<std::string, std::set<size_t> > &inc_dict, strvec &init_nodes, Algorithm& a){   
    std::cout<<"local-core(clique)\n";
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t N = init_nodes.size();
    intvec pcore(N); //
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    strInthashMap node_index; // (use hashtable instead of dictionary => Faster on large |V| datasets.)
    for(size_t i = 0; i<N; i++) node_index[init_nodes[i]] = i; // initialize node_index
    intvec nbrsizes(N); //
    // std::map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours.
    std::unordered_map<std::string, strset > init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
    std::vector<intvec> nbrs(N);
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(auto elem: e_id_to_edge){
        auto edge_sz = elem.second.size();
        sz_inc_edge += edge_sz;
        for(auto v_id: elem.second){
            auto j = node_index[v_id];  // j is of type int
            llb[j] = std::max(edge_sz - 1,llb[j]);
        // initialise number of neighbours and set of neighbours
            if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp.insert(u);
                        _tmp_sz+=1;
                    }
                }
                init_nbr[v_id] = std::move(_tmp);
                // init_nbrsize[v_id] =  _tmp.size();
                nbrsizes[j] = _tmp_sz;
                // sz_init_nbrs = _tmp_sz;
            }
            else{  // v_id exists in init_nbr map
                auto _tmp = &init_nbr[v_id];
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                // init_nbrsize[v_id] = init_nbr[v_id].size();
                // sz_init_nbrs -= nbrsizes[j];
                nbrsizes[j] = init_nbr[v_id].size();
                // sz_init_nbrs += nbrsizes[j];
            }
        }
    }
    for(auto i : init_nodes){
        sz_init_nbrs += init_nbr[i].size();
    }
    // std::cout<<"Init nbrs "<<sz_init_nbrs<<" Inc_edges "<<sz_inc_edge<<"\n"; 
    end3 = clock();
    // std::cout<<"Time for init_nbr calculation "<<double(end3 - start3) / double(CLOCKS_PER_SEC)<<"\n";
    time_t start4,end4;
    start4 = clock();

    // size_t nbrs_N[N+1];
    // nbrs_N[0] = 0;
    // size_t nbrs_F[sz_init_nbrs];
   
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    nbrs_N[0] = 0;
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));

    size_t glb = std::numeric_limits<size_t>::max();
    size_t gub = std::numeric_limits<size_t>::min();
    // #pragma omp parallel for default(none) shared(init_nodes,node_index,init_nbr,nbrs)
    int _i = 1, _index=0;
    for (auto node : init_nodes){
        nbrs_N[_i] = nbrs_N[_i-1] + init_nbr[node].size();
        _i++;
        size_t sz = 0;
        for (auto u: init_nbr[node]){
            // _tmp.push_back(node_index[u]);
            nbrs_F[_index++] = node_index[u];
            sz+=1;
        }
        // nbrs[node_index[node]] = std::move(_tmp);
        glb = std::min(glb, sz);
        gub = std::max(gub, sz);
    }
    // nbrs_N[N]++;
    end4 = clock();
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();
    // initialise core to a upper bound
    // #pragma omp parallel for default(none) shared(nbrsizes,pcore,N,glb,llb)
    for (size_t i = 0; i < N; i++){
        // printf("%d processing: %d\n",omp_get_thread_num(),i);
        pcore[i] = nbrsizes[i]; // initialize pcore
        // llb[i] = std::max(llb[i],glb);
    }
    
    if (log){
        strstrMap h0;
        for(size_t i=0;i<N;i++){
            h0[init_nodes[i]] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    // std::vector<size_t> hn(N);
    size_t iterations = 0;
    size_t correction_number=0, check=0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    // start_main = clock();
    while (1){
        iterations+=1;
        bool flag = true;
        // compute h-index and update core
        for(size_t i = 0; i<N; i++){
            // if (pcore[i] == llb[i]) continue;
            size_t H_value = hIndex_CSR(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if(H_value!=pcore[i])
            flag = false;
            pcore[i] = H_value;     //pcore[i] is same as hvn here
        }
       
        // end1 = clock();
        // a.core_exec_time += double(end1 - start1) / double(CLOCKS_PER_SEC);
        if (flag)
            break;

    }
    // end_main = clock();

    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        // a.nu_cu += nbrsizes[i] - pcore[i];
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);
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
            // std::cout<<"argv[2]"<<argv[2]<<"\n";
            if (argc>=3){
                if(std::string(argv[2])=="rand"){
                    std::cout<<"Generating random hyp.\n";
                    // getrandomHg(h);
                    getrandomHg(h,8,5,5);
                    h.dataset = "rand";
                }
                else{
                    getHg(argv[2],h);
                    h.dataset = argv[2];
                }
            }
            else{
                getHg("default",h);
                h.dataset = "default";
            }
            std :: string alg = "clique";
            if(argc>=4){
                alg = argv[3];
            }
            
            // std::cout<<argv[2]<<" "<<argv[3]<<"\n";
            std::cout<<"Clique graph\n";
            std::string init_type = "nbr"; // or "lub" (local upper bound)
            // h.initialise();
            clock_t ck_start = clock();
            getClique(h);
            h.initialise();
            if (argc>=5){
                if (atoi(argv[4])!=0)   h.writeneighborhood("../output/log_"+h.dataset+"_Nv.csv");
            }
            // h.printHypergraph();
            auto ck_time = double(clock()-ck_start)/double(CLOCKS_PER_SEC);
            // std::cout <<"E-Peel\n";
            Algorithm a(h);
            // if(alg == "E-Peel")
            // EPeel(h.dataset, h.e_id_to_edge, h.inc_dict, h.init_nodes, a);
            // else if(alg == "Local-core")
            local_core_clique(h.dataset, h.e_id_to_edge, h.inc_dict, h.init_nodes, a);
            a.output["algo"] = alg;
            a.output["num_threads"] = std::to_string(num_threads);
            a.output["total iteration"] = std::to_string(1);
            // a.output["execution time"] += ck_time;
            a.output["init_time"] += ck_time;
            a.write_results();
            std::cout<<"Execution time "<< a.exec_time<<"\n";
            // a.writecore();
            // a.writelog();
            // a.writeNbrQ();
            // std::cout << "core: \n";
            std::string file = "../output/core_"+a.output["algo"]+"_"+h.dataset+".csv";
            std::stringstream ss;
            for(auto elem: a.core)
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
    }
}
// g++ -std=c++11 -g -o ckmain clique_graph.cpp hypergraph.cpp  algorithms.cpp readhg.h