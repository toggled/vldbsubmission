from hgDecompose import utils
from time import time
import math
from copy import deepcopy



class HGDecompose():
    def __init__(self):
        self.bucket = {}
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


    def preprocess(self):
        pass

    def naiveNBR(self, H, verbose = True):
        start_execution_time = time()

        nodes = list(H.nodes)
        num_nodes = len(nodes)

        # Initial bucket fill-up
        for node in nodes:
            neighbors = list(H.neighbors(node))
            len_neighbors = len(neighbors)  # this computation can be repeated
            # node_to_neighbors[node] = neighbors
            self._node_to_num_neighbors[node] = len_neighbors
            # print(node, neighbors)
            if len_neighbors not in self.bucket:
                self.bucket[len_neighbors] = [node]
            else:
                self.bucket[len_neighbors].append(node)


        if(verbose):
            print("\n---------- Initial neighbors -------")
            for node in H.nodes:
                print(node, H.neighbors(node))
            print()

            print("\n---------- Initial bucket -------")
            print(self.bucket)
            print()

        

        for k in range(1, num_nodes + 1):

            while len(self.bucket.get(k, [])) != 0:
                v = self.bucket[k].pop(0)  # get first element in the

                if(verbose):
                    print("k:", k, "node:", v)
    
                self.core[v] = k
                temp_nodes = nodes[:] 
                temp_nodes.remove(v) 

                start_subgraph_time = time()
                H_temp = utils.strong_subgraph(H, temp_nodes) # Store.... + executation time.. 
                self.subgraph_time += time() - start_subgraph_time
                self.num_subgraph_call += 1

                # enumerating over all neighbors of v
                for u in utils.get_nbrs(H, v):  

                    if(verbose):
                        print(self._node_to_num_neighbors)
                        print("Considering neighbor", u)

                    start_neighborhood_call = time()
                    len_neighbors_u = utils.get_number_of_nbrs(H_temp, u)
                    self.neighborhood_call_time  += time() - start_neighborhood_call
                    self.num_neighborhood_computation += 1
                    # How many times is neighborhood computation done? and executation time...

                    max_value = max(len_neighbors_u, k)

                    if(verbose):
                        print("max core between", k, 'and', len_neighbors_u, "is ", max_value)
                        print("The location of", u, "is updated from", self._node_to_num_neighbors[u], "to", max_value)


                    # Move u to new location in bucket
                    start_bucket_update = time()
                    self.bucket[self._node_to_num_neighbors[u]].remove(u)
                    if max_value not in self.bucket:
                        # TODO does code reach here?
                        self.bucket[max_value] = [u]
                    else:
                        self.bucket[max_value].append(u)
                    self.num_bucket_update += 1
                    self.bucket_update_time += time() - start_bucket_update
                        
                    # How many times is bucket updated + executation time??? Store...

                    self._node_to_num_neighbors[u] = max_value

                    if(verbose):
                        print("-------- Updated bucket ---------")
                        print(self.bucket)
                        print()

                nodes = temp_nodes
                H = H_temp


        self.execution_time = time() - start_execution_time


    def naiveDeg(self, H, verbose = True):
        start_execution_time = time()

        nodes = list(H.nodes)
        num_nodes = len(nodes)

        # Initial bucket fill-up
        max_degree = -1
        for node in nodes:

            degree = H.degree(node)
            if(degree > max_degree):
                max_degree = degree
            self._node_to_degree[node] = degree
            # print(node, neighbors)
            if degree not in self.bucket:
                self.bucket[degree] = [node]
            else:
                self.bucket[degree].append(node)


        if(verbose):
            print("\n---------- Initial neighbors -------")
            for node in H.nodes:
                print(node, H.neighbors(node))
            print()

            print("\n---------- Initial bucket -------")
            print(self.bucket)
            print()

        

        for k in range(1, max_degree + 1):

            while len(self.bucket.get(k, [])) != 0:
                v = self.bucket[k].pop(0)  # get first element in the

                if(verbose):
                    print("k:", k, "node:", v)
    
                self.core[v] = k
                temp_nodes = nodes[:] 
                temp_nodes.remove(v) 

                start_subgraph_time = time()
                H_temp = utils.strong_subgraph(H, temp_nodes) # Store.... + executation time.. 
                self.subgraph_time += time() - start_subgraph_time
                self.num_subgraph_call += 1



                # enumerating over all neighbors of v
                for u in utils.get_nbrs(H, v):  

                    if(verbose):
                        print(self._node_to_degree)
                        print("Considering neighbor", u)

                    start_degree_call = time()
                    degree_u = utils.get_degree(H_temp, u)
                    self.degree_call_time  += time() - start_degree_call
                    self.num_degree_computation += 1
                    # How many times is neighborhood computation done? and executation time...

                    max_value = max(degree_u, k)

                    if(verbose):
                        print("max core between", k, 'and', degree_u, "is ", max_value)
                        print("The location of", u, "is updated from", self._node_to_num_neighbors[u], "to", max_value)


                    # Move u to new location in bucket
                    start_bucket_update = time()
                    self.bucket[self._node_to_degree[u]].remove(u)
                    if max_value not in self.bucket:
                        # TODO does code reach here?
                        self.bucket[max_value] = [u]
                    else:
                        self.bucket[max_value].append(u)
                    self.num_bucket_update += 1
                    self.bucket_update_time += time() - start_bucket_update
                        
                    # How many times is bucket updated + executation time??? Store...

                    self._node_to_degree[u] = max_value

                    if(verbose):
                        print("-------- Updated bucket ---------")
                        print(self.bucket)
                        print()

                nodes = temp_nodes
                H = H_temp


        self.execution_time = time() - start_execution_time

    # Interval generator function (s is a parameter)
    def generate_intervals(self, llb, lub, s = 1, verbose = False):
        min_llb = min([llb[u] for u in llb])
        ub_set = set([lub[u] for u in lub]).union([min_llb - 1])
        sorted_ub_set = sorted(ub_set, reverse=True)
        if(verbose):
            print('set of distinct values: ',sorted_ub_set)
        if s >= len(ub_set):
            yield sorted_ub_set[-1] + 1, sorted_ub_set[0]
        else:
            i = s
            while i < len(ub_set):
                yield sorted_ub_set[i] + 1, sorted_ub_set[i - s]
                if i+s < len(ub_set):
                    i += s
                else:
                    if i != len(ub_set) - 1:
                        yield sorted_ub_set[-1] + 1, sorted_ub_set[i]
                    i += s

    def improvedNBR(self, H, verbose = True):
        start_execution_time = time()
        nodes = list(H.nodes)
        num_nodes = len(nodes)
        

        
        lb1 = math.inf
        lb2 = {}
        ub1 = -math.inf
        # flag = {}
        # Initial bucket fill-up
        for node in nodes:
            nbrs = H.neighbors(node)
            # neighbors = list(H.neighbors(node))
            len_neighbors = len(nbrs) # this computation can be repeated
            lb1 = min(lb1,len_neighbors)
            ub1 = max(ub1,len_neighbors)
            # LB2 computation
            for u in nbrs:
                lb2[node] = min(lb2.get(node,len_neighbors),len(H.neighbors(u))-1)
            # node_to_neighbors[node] = neighbors
            self._node_to_num_neighbors[node] = len_neighbors
            # print(node, neighbors)
            if len_neighbors not in self.bucket:
                self.bucket[len_neighbors] = [node]
            else:
                self.bucket[len_neighbors].append(node)
            # flag[node] = False

        if(verbose):
            print("\n---------- Initial neighbors -------")
            for node in H.nodes:
                print(node, H.neighbors(node))
            print()


            print("\n---------- Initial bucket -------")
            print(self.bucket)
            print()

        for k in range(lb1, ub1):
            while len(self.bucket.get(k,[])) != 0:
                v = self.bucket[k].pop(0) # get first element in the
                if(verbose):
                    print("k:", k, "node:", v)
                self.core[v] = k

                temp_nodes = nodes[:] # Make a copy of nodes
                temp_nodes.remove(v)  # V' <- V \ {v}
                
                start_subgraph_time = time()
                H_temp = utils.strong_subgraph(H, temp_nodes) # Store.... + executation time.. 
                self.subgraph_time += time() - start_subgraph_time
                self.num_subgraph_call += 1
                
                
                for u in utils.get_nbrs(H, v):

                    if lb2[u] <= k:
                        start_neighborhood_call = time()
                        len_neighbors_u = utils.get_number_of_nbrs(H_temp, u)
                        self.neighborhood_call_time  += time() - start_neighborhood_call
                        self.num_neighborhood_computation += 1
                    

                        
                        max_value = max(len_neighbors_u, k)
                        if(verbose):
                            print("max core between", k, 'and', len_neighbors_u, "is ", max_value)
                            print("The location of", u, "is updated from", self._node_to_num_neighbors[u], "to", max_value)
                        
                        start_bucket_update = time()
                        self.bucket[self._node_to_num_neighbors[u]].remove(u)
                        if(max_value not in self.bucket):
                            self.bucket[max_value] = [u]
                        else:
                            self.bucket[max_value].append(u)
                        self.num_bucket_update += 1    
                        self.bucket_update_time += time() - start_bucket_update
                    

                        # update new location of u
                        self._node_to_num_neighbors[u] = max_value

                    if(verbose):
                        print("-------- Updated bucket ---------")
                        print(self.bucket)
                        print()
                nodes = temp_nodes
                H = H_temp

        self.execution_time = time() - start_execution_time

        if(verbose):
            print("\n\nOutput")
            print(self.core)


    def improved2NBR(self, H, s = 1, verbose = True):
        """ 
        :param H -> Hypergraph
        :param s -> Integer, algorithm parameter. 
        """
        start_execution_time = time()

        nodes = list(H.nodes)
        num_nodes = len(nodes)
        

        glb = math.inf
        llb = {}
        gub = -math.inf
        lub = {}
        # flag = {}
        # Initial bucket fill-up
        for node in nodes:
            nbrs = H.neighbors(node)
            # neighbors = list(H.neighbors(node))
            len_neighbors = len(nbrs)  # this computation can be repeated
            glb = min(glb, len_neighbors)
            gub = max(gub, len_neighbors)
            # LB2 computation
            for u in nbrs:
                llb[node] = min(llb.get(node, len_neighbors), len(H.neighbors(u)) - 1)
            # node_to_neighbors[node] = neighbors
            self._node_to_num_neighbors[node] = len_neighbors
            # print(node, neighbors)
            if len_neighbors not in self.bucket:
                self.bucket[len_neighbors] = [node]
            else:
                self.bucket[len_neighbors].append(node)
            # flag[node] = False

        if(verbose):
            print("\n---------- Initial neighbors -------")
            for node in H.nodes:
                print(node, H.neighbors(node))
            print()

            print("\n---------- Initial bucket -------")
            print(self.bucket)
            print()

        # Compute Local upper bounds
        copy_bucket = deepcopy(self.bucket)
        inv_bucket = deepcopy(self._node_to_num_neighbors)

        for k in range(glb, gub):
            while len(copy_bucket.get(k, [])) != 0:
                v = copy_bucket[k].pop(0)
                lub[v] = k
                for u in utils.get_nbrs(H, v):
                    if u not in lub:
                        max_value = max(inv_bucket[u] - 1, k)
                        copy_bucket[inv_bucket[u]].remove(u)
                        if(max_value not in copy_bucket):
                            copy_bucket[max_value] = [u]
                        else:
                            copy_bucket[max_value].append(u)
                        inv_bucket[u] = max_value

        if(verbose):
            print('local upper bound: ')
            print(sorted(lub.items()))
            print('local lower bound: ')
            print(sorted(llb.items()))


        gen = self.generate_intervals(llb, lub, s = s, verbose = verbose)
        final_bucket = {}
        setlb = {}
        inv_bucket = {}
        
        for lower, upper in gen:
            if(verbose):
                print("Inverval [%d,%d]"%(lower, upper))
            V_kmin = [u for u in nodes if lub[u] >= lower]
            for u in V_kmin:
                if u in self.core:
                    max_val = max(lower-1, llb[u], self.core[u])
                else:
                    max_val = max(lower-1, llb[u])
                final_bucket[max_val] = final_bucket.get(max_val, [])+[u]
                inv_bucket[u] = max_val
                setlb[u] = True


            start_subgraph_time = time()
            H_kmin = utils.strong_subgraph(H, V_kmin)
            self.subgraph_time += time() - start_subgraph_time
            self.num_subgraph_call += 1


            for k in range(lower-1, upper+1):
                while len(final_bucket.get(k, [])) != 0:
                    v = final_bucket[k].pop(0)
                    # core[v] = k
                    if setlb[v]:
                        start_neighborhood_call = time()
                        num_nbrs_v = utils.get_number_of_nbrs(H_kmin,v)
                        self.neighborhood_call_time  += time() - start_neighborhood_call
                        self.num_neighborhood_computation += 1
                    
                        final_bucket[num_nbrs_v] = final_bucket.get(num_nbrs_v,[])+[v]
                        inv_bucket[v] = num_nbrs_v
                        setlb[v] = False
                    else:
                        if k >= lower:
                            self.core[v] = k
                            setlb[v] = True
                        _temp_nodes = V_kmin[:]  # Make a copy of nodes
                        _temp_nodes.remove(v)  # V' <- V \ {v}
                        
                        start_subgraph_time = time()
                        _temp_Hkmin = utils.strong_subgraph(H_kmin, _temp_nodes)
                        self.subgraph_time += time() - start_subgraph_time
                        self.num_subgraph_call += 1
                        

                        for u in utils.get_nbrs(H_kmin, v):
                            start_neighborhood_call = time()
                            len_neighbors_u = utils.get_number_of_nbrs(_temp_Hkmin, u)
                            self.neighborhood_call_time  += time() - start_neighborhood_call
                            self.num_neighborhood_computation += 1
                    
                            max_value = max(len_neighbors_u, k)
                            start_bucket_update = time()
                            if max_value != inv_bucket[u]:
                                final_bucket[max_value] = final_bucket.get(max_value,[])+[u]
                                final_bucket[inv_bucket[u]].remove(u)
                                inv_bucket[u] = max_value
                                self.num_bucket_update += 1
                            self.bucket_update_time += time() - start_bucket_update
                        
                        V_kmin = _temp_nodes
                        H_kmin = _temp_Hkmin

        self.execution_time = time() - start_execution_time

        if(verbose):
            print("\n\nOutput")
            print(self.core)






