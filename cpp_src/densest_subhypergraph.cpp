#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include "hypergraph.h"
#include "readhg.h"

void write_subhg(std::string dataset, std::string alg, std::map<size_t, strvec > e_id_to_edge, std::set<size_t> subhg){

    std::string file = "../output/density/"+dataset+"_"+alg+".csv";
    std::cout<<file<<"\n";
    std::stringstream ss;
    for(auto e_id: subhg)
    {
        size_t sz = 0;
        for(auto u:e_id_to_edge[e_id]){
            if(sz!=e_id_to_edge[e_id].size()-1)
            ss << u <<",";
            else ss<<u;
            sz++;
        }
        ss<<"\n";
    }
    std::ofstream out(file.c_str());
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();

}

void write_result(std::string dataset, std::string alg, double density_nbr, double density_deg){

    std::string file = "../output/density/result.csv";
    std::cout<<file<<"\n";
    std::stringstream ss;
    ss <<dataset<<","<<alg<<","<<density_nbr<<","<<density_deg<<"\n";
    std::ofstream out(file.c_str(),std::ios::app);
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();

}

double nbr_density(std::set<size_t> &subhg, std::map<size_t, strvec > &e_id_to_edge){

    strIntMap nbrsizes;
    std::unordered_map<std::string, strset > init_nbr;
    size_t num_nodes = 0, sum_nbrs = 0;
    for(auto e_id : subhg){
        for(auto v_id: e_id_to_edge[e_id]){
            if ( init_nbr.find(v_id) == init_nbr.end() ) { 
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: e_id_to_edge[e_id]){
                    if (u!=v_id){
                        _tmp.insert(u);
                        _tmp_sz+=1;
                    }
                }
                init_nbr[v_id] = std::move(_tmp);
                nbrsizes[v_id] = _tmp_sz;
            }
            else{  
                auto _tmp = &init_nbr[v_id];
                for (auto u: e_id_to_edge[e_id]){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                nbrsizes[v_id] = init_nbr[v_id].size();
            }
        }
    }

    for(auto sz : nbrsizes){
        num_nodes+=1;
        sum_nbrs+=sz.second;
    }
    double density = double(sum_nbrs)/num_nodes;
    return density;

}

double deg_density(std::set<size_t> &subhg, std::map<size_t, strvec > &e_id_to_edge){

    std::map<std::string, std::set<size_t> > inc_dict;
    strvec init_nodes;
    size_t num_nodes = 0, sum_deg = 0;
    for(auto e_id:subhg){
        for(auto v_id:e_id_to_edge[e_id]){
            if ( inc_dict.find(v_id) == inc_dict.end() ) { // v not find (first time)
                inc_dict[v_id] = std::set<size_t>{e_id};
                init_nodes.push_back(v_id);
            }
            else 
            inc_dict[v_id].insert(e_id);
        }
    }

    for(auto u:init_nodes){
        num_nodes+=1;
        sum_deg+=inc_dict[u].size();
    }
    double density = double(sum_deg)/num_nodes;
    return density;


}

size_t get_number_of_nbrs(std::string v, std::map<std::string, std::set<size_t> > & inc_dict, std::map<size_t, strvec > &e_id_to_edge){
    std::set<std::string> nbrs;
    auto incident_edges = inc_dict[v];  //# {O(1)}
    if (incident_edges.size()){
        strboolMap visited_dict;
        for (auto e_id : incident_edges){  //# { O(d(v)) }
            for (auto u : e_id_to_edge[e_id]){  //# { O(|e|)}
                if (u != v){
                    if (visited_dict.find(u) == visited_dict.end()){ // u does not exists
                        visited_dict[u] = true;
                        nbrs.insert(u);
                    }
                    if (!visited_dict[u]){
                        visited_dict[u] = true;
                        nbrs.insert(u);
                    }
                }
            }
        }
    }
    return nbrs.size();
}

void removeV_transform(std::string v, std::map<std::string, std::set<size_t> > & inc_dict, std::map<size_t, strvec > &e_id_to_edge, std::set<std::string>& affected_nodes, std::set<size_t> &affected_edges){
    // For every edge e_id incident on v, remove e_id from every vertex u in edge[e_id] distinct from v
    for (auto e_id : inc_dict[v]){
        affected_edges.insert(e_id);
        for(auto u: e_id_to_edge[e_id]){
            if (u==v)   continue;
            inc_dict[u].erase(e_id);
            affected_nodes.insert(u);
        }
    }
    inc_dict[v] = std::set<size_t>();
}

void nbr_densest_hg(std::string dataset, std::map<size_t, strvec > e_id_to_edge, std::map<std::string, std::set<size_t> > inc_dict, strvec init_nodes, std::set<size_t> &subhg ){

    strIntMap nbrsizes;
    std::unordered_map<std::string, strset > init_nbr;
    std::map<int, std::set<std::string>> bucket;
    std::set<std::string> initnodes, affected_nodes;
    std::set<size_t> affected_edges; 
    size_t sum_nbrs = 0, num_nodes = init_nodes.size(), min_nbrs = init_nodes.size();
    double density = 0; 

    for(auto elem: e_id_to_edge){
        subhg.insert(elem.first);
        for(auto v_id: elem.second){
            if ( init_nbr.find(v_id) == init_nbr.end() ) { 
                auto _tmp = strset();
                int _tmp_sz = 0;
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp.insert(u);
                        _tmp_sz+=1;
                    }
                }
                init_nbr[v_id] = std::move(_tmp);
                nbrsizes[v_id] = _tmp_sz;
            }
            else{  
                auto _tmp = &init_nbr[v_id];
                for (auto u: elem.second){
                    if (u!=v_id){
                        _tmp->insert(u);
                    }
                }
                nbrsizes[v_id] = init_nbr[v_id].size();
            }
        }
    }

    for (auto node : init_nodes){
        initnodes.insert(node);
        auto len_neighbors = nbrsizes[node];
        if (bucket.find(len_neighbors) == bucket.end()){
            bucket[len_neighbors] = strset();
        }
        bucket[len_neighbors].insert(node);
        sum_nbrs += len_neighbors;
        min_nbrs = std::min(min_nbrs,len_neighbors);
    }

    density = double(sum_nbrs)/num_nodes;

    while(!initnodes.empty()){

        while(bucket[min_nbrs].size()==0)
        min_nbrs++;

        auto it = bucket[min_nbrs].begin();
        auto v = *it;
        // std::cout<<v<<" ";
        bucket[min_nbrs].erase(it);
        initnodes.erase(v);
        removeV_transform(v,inc_dict,e_id_to_edge,affected_nodes,affected_edges);

        for(auto u:affected_nodes){
            size_t len_nbrs = get_number_of_nbrs(u,inc_dict,e_id_to_edge);
            if(len_nbrs==0){

                bucket[nbrsizes[u]].erase(u);
                initnodes.erase(u);
                num_nodes--;
                sum_nbrs -= nbrsizes[u];
                nbrsizes[u] = 0;

            }else{

                bucket[nbrsizes[u]].erase(u);
                sum_nbrs -= nbrsizes[u];
                nbrsizes[u] = len_nbrs;
                bucket[nbrsizes[u]].insert(u);    
                sum_nbrs += nbrsizes[u];
                if(nbrsizes[u]<min_nbrs)
                min_nbrs = nbrsizes[u];            

            }

        }
        affected_nodes.clear();
        num_nodes--;
        sum_nbrs -= nbrsizes[v];
        // std::cout<<num_nodes<<" "<<sum_nbrs<<" ";
        if(num_nodes==0)
        break;
        double d2 = double(sum_nbrs)/num_nodes;
        // std::cout<<d2<<"\n";
        if(d2>density){
            density = d2;
            for(auto e_id:affected_edges){
                subhg.erase(e_id);
            }
            affected_edges.clear();
        }

    }

    std::cout<<"Maximum neighbourhood based density "<<density<<"\n";
    double density_deg = deg_density(subhg,e_id_to_edge);
    std::cout<<"Corresponding degree based density "<<density_deg<<"\n";
    write_result(dataset,"nbr",density,density_deg);

}

void deg_densest_hg(std::string dataset,std::map<size_t, strvec > e_id_to_edge, std::map<std::string, std::set<size_t> > inc_dict, strvec init_nodes, std::set<size_t> &subhg ){

    strIntMap degree;
    std::map<int, std::set<std::string>> bucket;
    std::set<std::string> initnodes, affected_nodes;
    std::set<size_t> affected_edges; 
    size_t sum_deg = 0, num_nodes = init_nodes.size(), min_deg = e_id_to_edge.size();
    double density = 0; 

    for(auto elem:e_id_to_edge)
    subhg.insert(elem.first);

    for (auto node : init_nodes){
        initnodes.insert(node);
        auto deg = inc_dict[node].size();
        degree[node] = deg;
        if (bucket.find(deg) == bucket.end()){
            bucket[deg] = strset();
        }
        bucket[deg].insert(node);
        sum_deg += deg;
        min_deg = std::min(min_deg,deg);
    }

    density = double(sum_deg)/num_nodes;

    while(!initnodes.empty()){

        while(bucket[min_deg].size()==0)
        min_deg++;

        auto it = bucket[min_deg].begin();
        auto v = *it;
        // std::cout<<v<<" ";
        bucket[min_deg].erase(it);
        initnodes.erase(v);
        removeV_transform(v,inc_dict,e_id_to_edge,affected_nodes,affected_edges);

        for(auto u:affected_nodes){
            size_t deg = inc_dict[u].size();
            if(deg==0){

                bucket[degree[u]].erase(u);
                initnodes.erase(u);
                num_nodes--;
                sum_deg -= degree[u];
                degree[u] = 0;

            }else{

                bucket[degree[u]].erase(u);
                sum_deg -= degree[u];
                degree[u] = deg;
                bucket[degree[u]].insert(u);    
                sum_deg += degree[u];
                if(degree[u]<min_deg)
                min_deg = degree[u];            

            }

        }
        affected_nodes.clear();
        num_nodes--;
        sum_deg -= degree[v];
        // std::cout<<num_nodes<<" "<<sum_deg<<" "<<subhg.size()<<" ";
        if(num_nodes==0){
            // std::cout<<"Hello\n";
            break;
        }
        double d2 = double(sum_deg)/num_nodes;
        // std::cout<<d2<<"\n";
        if(d2>density){
            density = d2;
            for(auto e_id:affected_edges){
                subhg.erase(e_id);
            }
            affected_edges.clear();
        }

    }
    // std::cout<<subhg.size()<<"\n";
    std::cout<<"Maximum degree based density "<<density<<"\n";
    double density_nbr = nbr_density(subhg,e_id_to_edge);
    std::cout<<"Corresponding neighbourhood based density "<<density_nbr<<"\n";
    write_result(dataset,"deg",density_nbr,density);

}

int main(int argc, char *argv[])
{
    
    std::string alg;
    Hypergraph h;
    if (argc>=2){
        getHg(argv[1],h);
        h.dataset = argv[1];
    }
    else{
        std::cout<<"No dataset provided. Exiting.";
        return 1;
    }
    if(argc>=3){
        alg = argv[2];
    }else{
        alg = "nbr";
    }
    h.initialise();
    std::set<size_t> subhg;
    if(alg=="nbr"){
        nbr_densest_hg(h.dataset,h.e_id_to_edge,h.inc_dict,h.init_nodes,subhg);
        write_subhg(h.dataset,alg,h.e_id_to_edge,subhg);
    } 
    if(alg=="deg"){
        deg_densest_hg(h.dataset,h.e_id_to_edge,h.inc_dict,h.init_nodes,subhg);
        write_subhg(h.dataset,alg,h.e_id_to_edge,subhg);
    }
    return 0;   
}

