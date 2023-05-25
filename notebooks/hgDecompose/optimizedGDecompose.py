from time import time
import math
from hgDecompose.Hypergraph import Hypergraph
from copy import deepcopy
from multiprocessing import Pool
from hgDecompose.utils import operator_H, par_operator_H
from hgDecompose.heapdict import heapdict
import pandas as pd
class GrDecompose():
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
    
    def generate_intervals(self, H, s = 1, verbose = False):
        sorted_ub_set = H.sorted_ub_set # precomputed
        len_ub_set = len(sorted_ub_set)
        if(verbose):
            print('set of distinct values: ',sorted_ub_set)
            print("#distinct values: ", len(sorted_ub_set))
            print('range: ',(sorted_ub_set[-1],sorted_ub_set[0]) )
        if s >= len_ub_set:
            yield sorted_ub_set[-1] + 1, sorted_ub_set[0]
        else:
            i = s
            while i < len_ub_set:
                yield sorted_ub_set[i] + 1, sorted_ub_set[i - s]
                if i+s < len_ub_set:
                    i += s
                else:
                    if i != len_ub_set - 1:
                        yield sorted_ub_set[-1] + 1, sorted_ub_set[i]
                    i += s
    def Graph_Core_decomp(self, H, V, lb1, ub1, setlb, bucket, inverse_bucket, verbose):
        if (verbose):
            print("Top-down Graph Core decomposition ")
            print("\n---------- Initial bucket -------")
            print(bucket)
            print()

            print('-- ',lb1,ub1,'---')
        num_nodes = 0
        _node_to_num_neighbors = {}
        bucket = {}
        for node in V:
            len_neighbors = H.get_number_of_nbrs(node)
            _node_to_num_neighbors[node] = len_neighbors
            if len_neighbors not in bucket:
                bucket[len_neighbors] = set()
            bucket[len_neighbors].add(node)
            num_nodes += 1

        if (verbose):
            print("\n Constructed Bucket: ")
            print(bucket)
            print()

        if (verbose):
            print('num_nodes: ',num_nodes)
        for k in range(1, num_nodes + 1):
            while len(bucket.get(k, [])) != 0:
                v = bucket[k].pop()  # get first element in the
                if(verbose):
                    print("k:", k, "node:", v)
            
                if k>=lb1:
                    self.core[v] = k
                nbr_v = H.neighbors(v)
                H.removeV_transform(v, False)

                # enumerating over all neighbors of v
                for u in nbr_v:
                    if (verbose):
                        print("node_to_num_neighbours: ",_node_to_num_neighbors)
                        print("Considering neighbor", u)

                    len_neighbors_u = H.get_number_of_nbrs(u)
                    max_value = max(len_neighbors_u, k)

                    if(verbose):
                        print("max core between", k, 'and', len_neighbors_u, "is ", max_value)
                        print("The location of", u, "is updated from", _node_to_num_neighbors[u], "to", max_value)


                    # Move u to new location in bucket
                    prev_idx = _node_to_num_neighbors[u]
                    bucket[prev_idx].remove(u)
                    if max_value not in bucket:
                        bucket[max_value] = set()
                    bucket[max_value].add(u)
                    _node_to_num_neighbors[u] = max_value

                    if(verbose):
                        print("-------- Updated bucket ---------")
                        print(bucket)
                        print()

    def Graph_top_down(self, H, s = 1, verbose = True):
        """ 
        Assume input H is a graph.
        This is not-so-efficient, yet correct implementation of a top-down algorithm.        
        """
        start_execution_time = time()
        llb = H.llb
        lub = H.lub

        start_init_time = time()
        Intervals = [] 
        gen = self.generate_intervals(H, s = s, verbose = verbose)
        print("*** ",H.llb," ***")
        for lower, upper in gen:
            Intervals.append((lower,upper))

        final_bucket = {}
        setlb = {}
        inv_bucket = {}
        
        start_loop_time = time()
        for lower,upper in Intervals:
            if(verbose):
                print("Inverval [%d,%d]"%(lower, upper))
            V_kmin = [u for u in H.init_node_iterator() if lub[u] >= lower]
            H_kmin = H.strong_subgraph(V_kmin)
            if(verbose):
                print("V_kmin: ", V_kmin)
                print("H[V_kmin]: ", H_kmin)

            self.Graph_Core_decomp(H_kmin, V_kmin, lower, upper, setlb, final_bucket, inv_bucket, verbose)
            if verbose:
                print('Partial Core: ',self.core)
                print("Bucket: ",final_bucket)
                print("setLB: ",setlb)
                print("=================")
        self.loop_time = time() - start_loop_time
        self.execution_time = time() - start_execution_time
        
        if(verbose):
            print("\n\nOutput")
            print(self.core)

    def Graph_Core_decomp2(self, H, lb1, ub1, setlb, bucket, inverse_bucket, verbose):
        """ INCORRECT: """
        if (verbose):
            print("Top-down Graph Core decomposition ")
            print("\n---------- Initial bucket -------")
            print(bucket)
            print()

            print('-- ',lb1,ub1,'---')
        # num_nodes = 0
        # _node_to_num_neighbors = {}
        # bucket = {}
        # for node in V:
        #     len_neighbors = H.get_number_of_nbrs(node)
        #     # if len_neighbors>lb1 and node not in self.core:
        #     #     len_neighbors = lb1 
        #     _node_to_num_neighbors[node] = len_neighbors
        #     if len_neighbors not in bucket:
        #         bucket[len_neighbors] = set()
        #     bucket[len_neighbors].add(node)
        #     num_nodes += 1

        # if (verbose):
        #     print("\n Constructed Bucket: ")
        #     print(bucket)
        #     print()

        # if (verbose):
        #     print('num_nodes: ',num_nodes)
        _node_to_num_neighbors = inverse_bucket
        for k in range(lb1, ub1 + 1):
            while len(bucket.get(k, [])) != 0:
                v = bucket[k].pop()  # get first element in the
                nbr_v = H.neighbors(v)
                len_nbr_v = len(nbr_v)
                if(verbose):
                    print("k:", k, "node:", v, ' ',' |N(v)| = ',len_nbr_v)
                if setlb[v] is True:
                    # if len_nbr_v>=lb1:
                    len_nbr_v = max(len_nbr_v,k)
                    if len_nbr_v not in bucket:
                        bucket[len_nbr_v] = []
                    if v not in bucket[len_nbr_v]:
                        bucket[len_nbr_v].append(v)

                    # update new location of u
                    inverse_bucket[v] = len_nbr_v
                    setlb[v] = False
                    if (verbose):
                        print('setlb False')
                        print(bucket)

                else:
                    if (verbose):
                        print('setlb True')
                        print('assign c(',v,') = ',k)
                    # if k>=lb1:
                    # if len_nbr_v>=k:
                    self.core[v] = k
                    setlb[v] = True 
                    H.removeV_transform(v, False)

                    # enumerating over all neighbors of v
                    for u in nbr_v:
                        if not setlb[u]:
                            if (verbose):
                                print("node_to_num_neighbours: ",_node_to_num_neighbors)
                                print("Considering neighbor", u)

                            len_neighbors_u = H.get_number_of_nbrs(u)
                            max_value = max(len_neighbors_u, k)

                            if(verbose):
                                print("max core between", k, 'and', len_neighbors_u, "is ", max_value)
                                print("The location of", u, "is updated from", _node_to_num_neighbors[u], "to", max_value)


                            # Move u to new location in bucket
                            prev_idx = _node_to_num_neighbors[u]
                            bucket[prev_idx].remove(u)
                            if max_value not in bucket:
                                bucket[max_value] = []
                            bucket[max_value].append(u)
                            _node_to_num_neighbors[u] = max_value

                            if(verbose):
                                print("-------- Updated bucket ---------")
                                print(bucket)
                                print()


    def Graph_top_down2(self, H, s = 1, verbose = True):
        """ 
        INCORRECT:
        Assume input H is a graph.
        A more efficient variant of Graph_top_down()      
        """
        start_execution_time = time()
        llb = H.llb
        lub = H.lub

        start_init_time = time()
        Intervals = [] 
        gen = self.generate_intervals(H, s = s, verbose = verbose)
        print("*** ",H.llb," ***")
        for lower, upper in gen:
            Intervals.append((lower,upper))

        final_bucket = {}
        setlb = {}
        inv_bucket = {}
        
        start_loop_time = time()
        for lower,upper in Intervals:
            if(verbose):
                print("Inverval [%d,%d]"%(lower, upper))
            V_kmin = [u for u in H.init_node_iterator() if lub[u] >= lower]
            H_kmin = H.strong_subgraph(V_kmin)
            if(verbose):
                print("V_kmin: ", V_kmin)
                print("H[V_kmin]: ", H_kmin)

            for u in V_kmin:
                if u in self.core:
                    # max_val = max(lower-1, min(llb[u], self.core[u]))
                    max_val = self.core[u]
                else:
                    max_val = lower

                if max_val not in final_bucket:
                    # final_bucket[max_val] = set()
                    final_bucket[max_val] = []
                if u not in final_bucket[max_val]:
                    final_bucket[max_val].append(u)
                inv_bucket[u] = max_val
                setlb[u] = True
                
            self.Graph_Core_decomp2(H_kmin, lower, upper, setlb, final_bucket, inv_bucket, verbose)
            if verbose:
                print('Partial Core: ',self.core)
                print("Bucket: ",final_bucket)
                print("setLB: ",setlb)
                print("=================")
        self.loop_time = time() - start_loop_time
        self.execution_time = time() - start_execution_time
        
        if(verbose):
            print("\n\nOutput")
            print(self.core)

class HDecompose():
    def __init__(self):
        # self.bucket = {}
        self.core = {}
        self._node_to_num_neighbors = {} # inverse bucket
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
    
    def generate_intervals(self, H, s = 1, verbose = False):
        sorted_ub_set = H.sorted_ub_set # precomputed
        len_ub_set = len(sorted_ub_set)
        if(verbose):
            print('set of distinct values: ',sorted_ub_set)
            print("#distinct values: ", len(sorted_ub_set))
            print('range: ',(sorted_ub_set[-1],sorted_ub_set[0]) )
        if s >= len_ub_set:
            yield sorted_ub_set[-1] + 1, sorted_ub_set[0]
        else:
            i = s
            while i < len_ub_set:
                yield sorted_ub_set[i] + 1, sorted_ub_set[i - s]
                if i+s < len_ub_set:
                    i += s
                else:
                    if i != len_ub_set - 1:
                        yield sorted_ub_set[-1] + 1, sorted_ub_set[i]
                    i += s

    def Core_decomp(self, H, V, lb1, ub1, setlb, bucket, inverse_bucket, verbose):
        """ Basic Peeling from [1, N] => Peel-algorithm"""
        if (verbose):
            print("Top-down Graph Core decomposition ")
            print("\n---------- Initial bucket -------")
            print(bucket)
            print()

            print('-- ',lb1,ub1,'---')
        num_nodes = 0
        _node_to_num_neighbors = {}
        bucket = {}
        for node in V:
            len_neighbors = H.get_number_of_nbrs(node)
            _node_to_num_neighbors[node] = len_neighbors
            if len_neighbors not in bucket:
                bucket[len_neighbors] = set()
            bucket[len_neighbors].add(node)
            num_nodes += 1

        if (verbose):
            print("\n Constructed Bucket: ")
            print(bucket)
            print()

        if (verbose):
            print('num_nodes: ',num_nodes)
        for k in range(1, num_nodes + 1):
            while len(bucket.get(k, [])) != 0:
                v = bucket[k].pop()  # get first element in the
                if(verbose):
                    print("k:", k, "node:", v)
                if k>=lb1:
                    self.core[v] = k
                nbr_v = H.neighbors(v)
                H.removeV_transform(v, False)

                # enumerating over all neighbors of v
                for u in nbr_v:
                    if (verbose):
                        print("node_to_num_neighbours: ",_node_to_num_neighbors)
                        print("Considering neighbor", u)

                    len_neighbors_u = H.get_number_of_nbrs(u)
                    max_value = max(len_neighbors_u, k)

                    if(verbose):
                        print("max core between", k, 'and', len_neighbors_u, "is ", max_value)
                        print("The location of", u, "is updated from", _node_to_num_neighbors[u], "to", max_value)


                    # Move u to new location in bucket
                    prev_idx = _node_to_num_neighbors[u]
                    bucket[prev_idx].remove(u)
                    if max_value not in bucket:
                        bucket[max_value] = set()
                    bucket[max_value].add(u)
                    _node_to_num_neighbors[u] = max_value

                    if(verbose):
                        print("-------- Updated bucket ---------")
                        print(bucket)
                        print()

    def naive_top_down(self, H, s = 1, verbose = True):
        """ 
        This is not-so-efficient, yet correct implementation of a top-down algorithm.        
        """
        llb = H.llb
        lub = H.lub

        Intervals = [] 
        gen = self.generate_intervals(H, s = s, verbose = verbose)
        if (verbose):
            print("*** ",H.llb," ***")
        for lower, upper in gen:
            Intervals.append((lower,upper))

        final_bucket = {}
        setlb = {}
        inv_bucket = {}

        for lower,upper in Intervals:
            if(verbose):
                print("Inverval [%d,%d]"%(lower, upper))
            V_kmin = [u for u in H.init_node_iterator() if lub[u] >= lower]
            H_kmin = H.strong_subgraph(V_kmin)
            if(verbose):
                print("V_kmin: ", V_kmin)
                print("H[V_kmin]: ", H_kmin)

            self.Core_decomp(H_kmin, V_kmin, lower, upper, setlb, final_bucket, inv_bucket, verbose)
            if verbose:
                print('Partial Core: ',self.core)
                print("Bucket: ",final_bucket)
                print("setLB: ",setlb)
                print("=================")
        
        if(verbose):
            print("\n\nOutput")
            print(self.core)

    def Core_decomp2(self, H, V, lb1, ub1, setlb, bucket, inverse_bucket, verbose):
        """ 
            Better Peeling. 
            1. We peel only upto a point needed [1,ub1].
            2. If core number for a vertex u is already computed:
                2.1) initialize u to its core-number in the bucket during initialisation phase.
                2.2) After v \in Nbr(u) is assigned core-number, while we update nbr(v)'s bucket index, we can ignore v from such update. 
         """
        if (verbose):
            print("Top-down Graph Core decomposition ")
            print("\n---------- Initial bucket -------")
            print(bucket)
            print()

            print('-- ',lb1,ub1,'---')
        num_nodes = 0
        _node_to_num_neighbors = {}
        bucket = {}
        already_computed_flag = {}
        for node in V:
            if node in self.core:
                len_neighbors = self.core[node]
                already_computed_flag[node] = True 
            else:
                already_computed_flag[node] = False 
                len_neighbors = H.get_number_of_nbrs(node)
            _node_to_num_neighbors[node] = len_neighbors
            if len_neighbors not in bucket:
                bucket[len_neighbors] = set()
            bucket[len_neighbors].add(node)
            num_nodes += 1

        if (verbose):
            print("\n Constructed Bucket: ")
            print(bucket)
            print()

        if (verbose):
            print('num_nodes: ',num_nodes)
        for k in range(1, ub1+1):
            while len(bucket.get(k, [])) != 0:
                v = bucket[k].pop()  # get first element in the
                if(verbose):
                    print("k:", k, "node:", v)
                if k>=lb1:
                    self.core[v] = k
                    already_computed_flag[v] = True 
                nbr_v = H.neighbors(v)
                H.removeV_transform(v, False)

                # enumerating over all neighbors of v
                for u in nbr_v:
                    if not already_computed_flag[u]:
                        if (verbose):
                            print("node_to_num_neighbours: ",_node_to_num_neighbors)
                            print("Considering neighbor", u)

                        len_neighbors_u = H.get_number_of_nbrs(u)
                        max_value = max(len_neighbors_u, k)

                        if(verbose):
                            print("max core between", k, 'and', len_neighbors_u, "is ", max_value)
                            print("The location of", u, "is updated from", _node_to_num_neighbors[u], "to", max_value)


                        # Move u to new location in bucket
                        prev_idx = _node_to_num_neighbors[u]
                        bucket[prev_idx].remove(u)
                        if max_value not in bucket:
                            bucket[max_value] = set()
                        bucket[max_value].add(u)
                        _node_to_num_neighbors[u] = max_value

                        if(verbose):
                            print("-------- Updated bucket ---------")
                            print(bucket)
                            print()

    def top_down(self, H, s = 1, verbose = True):
        """ 
        This is not-so-efficient, yet correct implementation of a top-down algorithm.        
        """
        start_execution_time = time()
        llb = H.llb
        lub = H.lub

        Intervals = [] 
        gen = self.generate_intervals(H, s = s, verbose = verbose)
        if (verbose):
            print("*** ",H.llb," ***")
        for lower, upper in gen:
            Intervals.append((lower,upper))

        final_bucket = {}
        setlb = {}
        inv_bucket = {}

        for lower,upper in Intervals:
            if(verbose):
                print("Inverval [%d,%d]"%(lower, upper))
            V_kmin = [u for u in H.init_node_iterator() if lub[u] >= lower]
            H_kmin = H.strong_subgraph(V_kmin)
            if(verbose):
                print("V_kmin: ", V_kmin)
                print("H[V_kmin]: ", H_kmin)

            self.Core_decomp2(H_kmin, V_kmin, lower, upper, setlb, final_bucket, inv_bucket, verbose)
            if verbose:
                print('Partial Core: ',self.core)
                print("Bucket: ",final_bucket)
                print("setLB: ",setlb)
                print("=================")
        
        if(verbose):
            print("\n\nOutput")
            print(self.core)
        self.execution_time = time() - start_execution_time