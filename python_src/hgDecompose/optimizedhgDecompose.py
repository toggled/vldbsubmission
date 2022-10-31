from json.encoder import INFINITY
from time import time
import math
from hgDecompose.Hypergraph import Hypergraph
from copy import deepcopy
from multiprocessing import Pool
from hgDecompose.utils import operator_H, par_operator_H
from hgDecompose.heapdict import heapdict
import pandas as pd
# from tests.verify_kcore import *

class HGDecompose():
    def __init__(self):
        # self.bucket = {}
        self.core = {}
        # self._node_to_num_neighbors = {} # inverse bucket
        self._node_to_degree = {}
        self.execution_time = 0
        self.bucket_update_time = 0
        self.neighborhood_call_time = 0
        self.degree_call_time = 0
        self.num_neighborhood_computation = 0
        self.num_degree_computation = 0
        self.num_bucket_update = 0
        self.subgraph_time = 0
        self.num_subgraph_call = 0
        self.init_time = 0
        self.loop_time = 0
        self.inner_iteration = 0
        self.total_iteration = 0
        self.core_correct_time = 0
        self.h_index_time = 0
        self.max_n = 0 # For #iterations vs dataset barplot
        self.core_correctionvol_n = [] #  core_corrections volume per iteration => Ammount of core_correction done. => Relation with runtime
        self.core_correction_volume = 0 # For core_correction volume vs dataset plot
        self.reduction_hhat_n = [] # [ hhat^{n-1} - hhat^{n}, for n \in [1, tau] ] => Convergence plot.

    def preprocess(self):
        pass

    def naiveNBR(self, H, verbose = True):
        start_execution_time = time()
        num_nodes = 0
        _node_to_num_neighbors = {}
        bucket = {}
        # Initial bucket fill-up
        start_init_time = time()
        for node in H.init_node_iterator():
            len_neighbors = H.get_init_nbrlen(node)
            _node_to_num_neighbors[node] = len_neighbors
            if len_neighbors not in bucket:
                bucket[len_neighbors] = set()
            bucket[len_neighbors].add(node)
            num_nodes += 1
        self.init_time = time() - start_init_time

        if(verbose):
            # print("\n---------- Initial neighbors -------")
            # for node in H.nodes():
            #     print(node, H.neighbors(node))
            # print()

            print("\n---------- Initial bucket -------")
            print(bucket)
            print()

        
        start_loop_time = time()
        for k in range(1, num_nodes + 1):
            while len(bucket.get(k, [])) != 0:
                # print(k)
                v = bucket[k].pop()  # get first element in the
                
                if(verbose):
                    print("k:", k, "node:", v)
    
                self.core[v] = k
                
                start_neighborhood_call = time()
                nbr_v = H.neighbors(v)
                self.neighborhood_call_time += time() - start_neighborhood_call
                self.num_neighborhood_computation += 1

                start_subgraph_time = time()
                H.removeV_transform(v, False)
                # H.removeV_transform2(v,verbose)
                self.subgraph_time += time() - start_subgraph_time
                self.num_subgraph_call += 1

                # enumerating over all neighbors of v
                for u in nbr_v:
                    self.inner_iteration += 1
                    self.total_iteration +=1
                    if (verbose):
                        print("node_to_num_neighbours: ",_node_to_num_neighbors)
                        print("Considering neighbor", u)

                    start_neighborhood_call = time()
                    len_neighbors_u = H.get_number_of_nbrs(u)
                    self.neighborhood_call_time  += time() - start_neighborhood_call
                    self.num_neighborhood_computation += 1

                    max_value = max(len_neighbors_u, k)

                    if(verbose):
                        print("max core between", k, 'and', len_neighbors_u, "is ", max_value)
                        print("The location of", u, "is updated from", _node_to_num_neighbors[u], "to", max_value)


                    # Move u to new location in bucket
                    start_bucket_update = time()
                    prev_idx = _node_to_num_neighbors[u]

                    bucket[prev_idx].remove(u)

                    if max_value not in bucket:
                        bucket[max_value] = set()
                    bucket[max_value].add(u)
                    self.num_bucket_update += 1
                    self.bucket_update_time += time() - start_bucket_update
                    _node_to_num_neighbors[u] = max_value

                    if(verbose):
                        print("-------- Updated bucket ---------")
                        print(bucket)
                        print()
        # print(self.core)
        self.loop_time = time() - start_loop_time
        self.execution_time = time() - start_execution_time
        if (verbose):
            print("\n\nOutput")
            print(self.core)
            
    def naiveDeg(self, H, verbose = True):
        assert isinstance(H,Hypergraph)
        start_execution_time = time()
        bucket = {}
        nodes = H.init_nodes
        num_nodes = len(nodes)

        # Initial bucket fill-up
        max_degree = -1
        for node in nodes:

            degree = H.degree(node)
            if(degree > max_degree):
                max_degree = degree
            self._node_to_degree[node] = degree
            # print(node, neighbors)
            if degree not in bucket:
                bucket[degree] = [node]
            else:
                bucket[degree].append(node)


        if(verbose):
            print("\n---------- Initial neighbors -------")
            for node in H.nodes():
                print(node, H.neighbors(node))
            print()

            print("\n---------- Initial bucket -------")
            print(bucket)
            print()

        

        for k in range(1, max_degree + 1):
            while len(bucket.get(k, [])) != 0:
                v = bucket[k].pop(0)  # get first element in the

                if(verbose):
                    print("k:", k, "node:", v)
    
                self.core[v] = k
                # temp_nodes = nodes[:] 
                # temp_nodes.remove(v) 
                nbr_v = H.neighbors(v)

                start_subgraph_time = time()
                # H_temp = H.strong_subgraph(temp_nodes) # Store.... + executation time.. 
                H.removeV_transform(v, False)
                self.subgraph_time += time() - start_subgraph_time
                self.num_subgraph_call += 1



                # enumerating over all neighbors of v
                for u in nbr_v:  
                    if(verbose):
                        print(self._node_to_degree)
                        print("Considering neighbor", u)

                    start_degree_call = time()
                    degree_u = H.degree(u)
                    self.degree_call_time  += time() - start_degree_call
                    self.num_degree_computation += 1
                    # How many times is neighborhood computation done? and executation time...

                    max_value = max(degree_u, k)

                    if(verbose):
                        print("max core between", k, 'and', degree_u, "is ", max_value)
                        print("The location of", u, "is updated from", self._node_to_degree[u], "to", max_value)


                    # Move u to new location in bucket
                    start_bucket_update = time()
                    bucket[self._node_to_degree[u]].remove(u)
                    if max_value not in bucket:
                        # TODO does code reach here?
                        bucket[max_value] = [u]
                    else:
                        bucket[max_value].append(u)
                    self.num_bucket_update += 1
                    self.bucket_update_time += time() - start_bucket_update
                        
                    # How many times is bucket updated + executation time??? Store...

                    self._node_to_degree[u] = max_value

                    if(verbose):
                        print("-------- Updated bucket ---------")
                        print(bucket)
                        print()

                # nodes = temp_nodes
                # H = H_temp


        self.execution_time = time() - start_execution_time