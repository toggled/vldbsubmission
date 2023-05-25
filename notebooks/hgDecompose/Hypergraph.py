import math
import itertools
import random
from typing import final 

class Hypergraph:
    """ 
    Our own hypergraph representation class. 
    We store hyperedge list in compressed format using two things- 1) e_indices (a dict) 2) e_nodes (a list)
    Although edge-centric queries (e.g. edge enumeration) are facilitated in this way, node-centric queries are not convenient.
    To support node-centric queries, we also maintain incidence dictionary inc_dict (key = v_ids, values = incident edge ids)
    """

    def __init__(self, _edgedict=None):
        
        self.e_indices = {}  # (position, position+edge_size) of edge e in e_nodes list
        self.e_nodes = []  # flattened edge list
        self.inc_dict = {}  # key: nodeid, value = ids of incident edges (set)
        # degree pre-compute => degree_dict or len_incedge = {}
        self.degree_dict = {}
        self.init_nbrsize = {} # initial nbrhood sizes. can be precomputed.
        self.init_nbr = {}
        self.init_eids = {}
        self.init_nodes = []
        if _edgedict is None or len(_edgedict)==0:  # Returns an empty Hypergraph
            return

        self.i = 0
        for e_id, e in _edgedict.items():
            _len = len(e)
            
            self.e_indices[e_id] = (self.i, self.i + _len)
            self.init_eids[e_id] = (self.i, self.i + _len)
            for v in e:
                self.e_nodes.append(v)
                if v not in self.inc_dict:
                    self.inc_dict[v] = set()  # create incident edge entry for v
                    self.init_nodes.append(v)
                self.inc_dict[v].add(e_id)  # incident edge update
                self.degree_dict[v] = self.degree_dict.get(v, 0) + 1  # degree update
                nbr_v = self.init_nbr.get(v, set()).union(e)
                nbr_v.remove(v)
                self.init_nbrsize[v] = len(nbr_v)  # neighbourhood length update
                self.init_nbr[v] = nbr_v  # neighbourbood set update
            self.i += _len

        self.init_nodes = sorted(self.init_nodes)
        
        # self.compute_Bounds()
        self.compute_Bounds2()

    def compute_Bounds2(self):
        """ Computes Local Upper Bound differently from compute_Bounds() """
        # print('Global-local ub')
        # Computing global upper and lower bounds
        self.glb = math.inf
        self.gub = -math.inf
        # Computing local lower bounds
        # self.precomputedlb2 = {}
        # # Computing local upper bounds
        self.lub = {}
        
        self.sorted_ub_set = set() # || => upper bound for param_s
        for v in self.node_iterator():
            len_neighbors_v = self.init_nbrsize[v]
            self.glb = min(self.glb,len_neighbors_v)
            self.gub = max(self.gub, len_neighbors_v)
            self.lub[v] = len_neighbors_v
            self.sorted_ub_set.add(len_neighbors_v)

        # ## print('---- ',self.glb, self.gub)
        # for k in range(self.glb, self.gub+1):
        #     while len(_bucket.get(k, [])) != 0:
        #         v = _bucket[k].pop()
        #         self.lub[v] = k
        #         self.sorted_ub_set.add(k) # add computed local upper bound to ub_set
        #         for u in self.init_nbr[v]:
        #             if u not in self.lub:
        #                 max_value = max(_inv_bucket[u] - 1, k)
        #                 _bucket[_inv_bucket[u]].remove(u)
        #                 if(max_value not in _bucket):
        #                     _bucket[max_value] = set()
        #                 _bucket[max_value].add(u)
        #                 _inv_bucket[u] = max_value

        
        self.llb = {}
        _min_llb = math.inf
        # Local lower bound computation
        for v in self.node_iterator():
            _max = -math.inf
            for e_id in self.inc_dict[v]:
                _max = max(_max, len(self.get_edge_byindex(e_id)) - 1)
            self.llb[v] = max(_max, self.glb)
            _min_llb = min(_min_llb, self.llb[v])

        self.sorted_ub_set.add(_min_llb - 1)
        self.sorted_ub_set = sorted(list(self.sorted_ub_set), reverse=True)

    def compute_Bounds(self):
        # print('Global-local ub')
        # Computing global upper and lower bounds
        self.glb = math.inf
        self.gub = -math.inf
        # Computing local lower bounds
        # self.precomputedlb2 = {}
        # # Computing local upper bounds
        self.lub = {}
        
        # Auxiliary variables to assist computation of ub2
        _inv_bucket = {}
        _bucket = {}
        self.sorted_ub_set = set() # || => upper bound for param_s
        for v in self.node_iterator():
            len_neighbors_v = self.init_nbrsize[v]
            self.glb = min(self.glb,len_neighbors_v)
            self.gub = max(self.gub, len_neighbors_v)
            _inv_bucket[v] = len_neighbors_v
            ## print(node, neighbors)
            if len_neighbors_v not in _bucket:
                _bucket[len_neighbors_v] = set()
            _bucket[len_neighbors_v].add(v)
        ## print('---- ',self.glb, self.gub)
        for k in range(self.glb, self.gub+1):
            while len(_bucket.get(k, [])) != 0:
                v = _bucket[k].pop()
                self.lub[v] = k
                self.sorted_ub_set.add(k) # add computed local upper bound to ub_set
                for u in self.init_nbr[v]:
                    if u not in self.lub:
                        max_value = max(_inv_bucket[u] - 1, k)
                        _bucket[_inv_bucket[u]].remove(u)
                        if(max_value not in _bucket):
                            _bucket[max_value] = set()
                        _bucket[max_value].add(u)
                        _inv_bucket[u] = max_value

        
        self.llb = {}
        _min_llb = math.inf
        # Local lower bound computation
        for v in self.node_iterator():
            _max = -math.inf
            for e_id in self.inc_dict[v]:
                _max = max(_max, len(self.get_edge_byindex(e_id)) - 1)
            self.llb[v] = max(_max, self.glb)
            _min_llb = min(_min_llb, self.llb[v])

        self.sorted_ub_set.add(_min_llb - 1)
        self.sorted_ub_set = sorted(list(self.sorted_ub_set), reverse=True)
        del _bucket
        del _inv_bucket

    def get_init_nbr(self, v):
        return self.init_nbr[v]

    def get_init_nbrlen(self, v):
        return self.init_nbrsize[v]

    def add_edge(self, e_id, e_nodes):
        """ Add an edge to the hypergraph. It does not check repeated edge."""
        _len = len(e_nodes)
        self.e_indices[e_id] = (self.i, self.i + _len)
        for v in e_nodes:
            self.e_nodes.append(v)
            if v not in self.inc_dict:
                self.inc_dict[v] = set()  # create incident edge entry for v
                self.init_nodes.append(v)
            self.inc_dict[v].add(e_id)  # incident edge update
            self.degree_dict[v] = self.degree_dict.get(v, 0) + 1  # degree update

        self.i += _len

    def get_stronglyinduced_edgeIds(self, vertex_list):
        """ Returns the list of edge_ids strongly induced by vertex_list"""
        temp_H = self.strong_subgraph(vertex_list)
        assert isinstance(temp_H, Hypergraph)
        return [e_id for e_id in temp_H.e_indices.keys()]

    def del_edge(self, e_id):
        """ Delete an edge given a set of nodes e."""
        if e_id not in self.e_indices: # If edge e_id does not exist in the hypergraph, do nothing.
            return 
        for v in self.get_edge_byindex(e_id):
            self.inc_dict[v].remove(e_id) # Remove e_id from v's incident edge list.
            if len(self.inc_dict[v])==0: # After deleting e_id there is no edge incident on v anymore. Then del v
                del self.inc_dict[v]
            self.degree_dict[v] -= 1 # decrease degree of v by 1
        
        del self.e_indices[e_id]

    def hasEdge(self, e_nodes):
        Exists = False 
        for v in e_nodes:
            if v in self.inc_dict:
                for e_id in self.inc_dict[v]:
                    if tuple(self.get_edge_byindex(e_id)) == e_nodes:
                        Exists = True 
                        break 
            if Exists:
                break 
        return Exists 

    def sample_v_preferential_attachment(self, num_sample):
        """ 
        Sample num_sample vertices following degree distribution. 
        Algorithm: Reservoir sampling  (WeightedReservoir-Chao wikipedia: https://en.wikipedia.org/wiki/Reservoir_sampling)
        """
        # WeightedReservoir-Chao (S[1..n], R[1..k])
        WSum = 0
        R = []
        #fill the reservoir array
        k = 0
        n = len(self.init_nodes)
        # print(self.init_nodes,' num_sample: ',num_sample)
        for i in range(num_sample):
            R.append(self.init_nodes[i])
            WSum = WSum + self.degree(self.init_nodes[i])
            k+=1
        
        # print(R)
        for i in range(k, n):
            # print(i,'-', R)
            WSum = WSum + self.degree(self.init_nodes[i])
            p = self.degree(self.init_nodes[i])*1.0 / WSum # probability for this item
            j = random.random();          # uniformly random between 0 and 1
            if j <= p:               # select item according to probability
                R[random.randint(0,k-1)] = self.init_nodes[i]  #uniform selection in reservoir for replacement
        return R 

    def get_edge_byindex(self, e_id):
        """ Return edge by edge_id """
        e_start, e_end = self.e_indices[e_id]
        return self.e_nodes[e_start:e_end]

    def edge_iterator(self):
        """ returns: iterator """
        for e_id in self.e_indices.keys():
            yield self.get_edge_byindex(e_id)

    def edge_eid_iterator(self):
        """ returns: iterator """
        for e_id in self.e_indices.keys():
            yield (e_id, self.get_edge_byindex(e_id))

    def init_node_iterator(self):
        """ 
        Returns: iterator of initial nodes. 
        Faster than node_iterator() and deterministic.
        """
        for v in self.init_nodes:
            yield v

    def node_iterator(self):
        """ returns: iterator """
        for v_id in self.inc_dict.keys():
            yield v_id

    def nodes(self):
        """ returns: list of vertices """
        return [v for v in self.node_iterator()]

    def edges(self):
        """ returns: list of edges (each edge is a list of vertex ids) """
        return [e for e in self.edge_iterator()]

    def degree(self, u):
        """ returns: integer """
        # assert (len(self.inc_dict.get(u,[])) == self.degree_dict[u])
        return self.degree_dict.get(u, 0)

    def dim(self, e):
        """ returns: integer """
        return len(e) - 1

    def neighbors(self, v):
        return [u for u in self.neighbors_iterator(v)]

    def get_number_of_nbrs(self, u):
        return len(self.neighbors(u))

    def compute_degree(self, u):
        return len(self.inc_dict[u])
        
    def degree_density(self):
        # sum_deg = 0
        # for u in self.inc_dict:
        #     sum_deg += self.compute_degree(u)
        # return sum_deg/self.get_N()
        return self.get_M()*1.0 / self.get_N()

    def neighbors_iterator(self, v):
        """ Returns the set of neighbours of v.
            implements a traversal from vertex v to each of its neighbours in contrast to set in neighbors(). 
            It also returns an iterator. So it avoids creating the neighborhood list explicitely.
            Overall complexity: O(d(v) * |e_max|), where e_max = largest hyperedge 
        """
        incident_edges = self.inc_dict.get(v, None)  # {O(1)}
        if incident_edges:
            visited_dict = {}
            for e_id in incident_edges:  # { O(d(v)) }
                for u in self.get_edge_byindex(e_id):  # { O(|e|)}
                    if u != v:
                        if not visited_dict.get(u, False):
                            visited_dict[u] = True
                            yield u
        else:
            return
            
    def strong_subgraph(self, vertex_list):
        """ returns: Hypergraph object. """
        H = Hypergraph()
        e_indices = {}  # (position, edge_size) of edge e in e_nodes list
        e_nodes = []  # flattened edge list
        inc_dict = {}
        H.i = 0
    
        # print('inc_dict: ',self.inc_dict.items())
        # print('e_indices: ',self.e_indices.items())
        # print('e_nodes: ',self.e_nodes)
    
        for e_id in self.e_indices:
            e = self.get_edge_byindex(e_id)
            flag = True 
            for u in e:
                if u not in vertex_list:
                    flag = False
                    break
            if flag:
                e_nodes += e
                _lene = len(e)
                e_indices[e_id] = (H.i, H.i + _lene)
                for v in e:
                    if v not in inc_dict:
                        inc_dict[v] = set()
                    if e_id not in inc_dict[v]:
                        inc_dict[v].add(e_id)
                        H.degree_dict[v] = H.degree_dict.get(v,0) + 1
                H.i += _lene
            # else:
            #     inc_dict[v] = inc_dict.get(v, []) + [e_id]
    
        # print('After: ','inc_dict = ',inc_dict.items(),'\n','e_indices = ',e_indices,'\n',' e_nodes = ',e_nodes)
        H.e_nodes = e_nodes
        H.inc_dict = inc_dict
        H.e_indices = e_indices
        return H

    def removeV_transform(self, v, verbose=False):
        """ removes input vertex v and transforms this hypergraph into a sub-hypergraph strongly induced by V\{v}
        Here we do not maintain nbr and len_nbr dictionaries.
        """
        incident_eids = set()  # set of edge_ids incident on v
        for e_id in self.inc_dict.get(v, []):
            incident_eids.add(e_id)

        if verbose:
            print("incident edges on ",v," : ", incident_eids)

        # Update incident edges and degree of every nbr of v
        for u in self.neighbors_iterator(v): # traverse over neighbours of v
            if verbose:
                print('traversing nbr: ',u)
            self.inc_dict[u] -= incident_eids # remove v's incident edges from u's incident edges.
            self.degree_dict[u] = len(self.inc_dict.get(u, []))

        if v in self.inc_dict:
            del self.inc_dict[v]

        if v in self.degree_dict:
            del self.degree_dict[v]
    # TODO
    # def addV_transform(self, S):
    #     pass
    #     prev_v = prev_v.union(S)
    #     for every edge e:
    #         if e \subset prev_v:
    #             if e.id not already exist:
    #                 add (e)

    def adV_transform(self, S):
        """ 
        S: is a set of vertices
        We assume the current hypergraph is already a strong subgraph. meaning the inc_dict and e_indices are maintained.
        """
        self.prev_V = self.prev_V.union(S)
        for e in self.init_eids: # But this will not give me all the edge_id's
            start_e,end_e = self.init_eids[e]
            edge = self.e_nodes[start_e:end_e]
            if set(edge).issubset(self.prev_V):
                self.e_indices[e] = self.init_eids[e]
                for u in edge:
                    self.inc_dict[u].add(e)
                
                
    # def removeV_transform(self, v, verbose=False):
    #     """ removes input vertex v and transforms this hypergraph into a sub-hypergraph strongly induced by V\{v} """
    #
    #     # inc_on_v = self.inc_dict.get(v,[]) # List of incident edges(v) { o(1) }
    #     # if len(inc_on_v): # If 0 incident edges, the hypergraph won't need any update.
    #     # for e_id in inc_on_v: # For every edges incident on v { O(d(v)) }
    #     # spos_e, epos_e = self.e_indices[e_id] # { o(1) }
    #     # del self.e_indices[e_id] # Delete that edge (to be precise, we remove its index only) { o(1) }
    #     # # Check if the neighbours of v in e needs to be removed as well.
    #     # for nbrv_in_e in self.e_nodes[spos_e:epos_e]: # For every u in e { O(|e_max|) < O(N) }
    #     #     if nbrv_in_e in self.inc_dict:
    #     #         if len(self.inc_dict[nbrv_in_e]) <= 1: # If u has no edge except e incident on it
    #     #             del self.inc_dict[nbrv_in_e] # Remove u
    #     #         else:
    #     #             self.inc_dict[nbrv_in_e].remove(e_id) # Otherwise, update the incidence entry of u
    #
    #     incident_edges = self.inc_dict.get(v, None)  # {O(1)}
    #     if (verbose):
    #         print(incident_edges)
    #     if incident_edges:
    #         visited_dict = {}
    #         # isolated_vertices = set([v]) # the set of nbrs of v who became isolated after removal of v. They can not be in any edge.
    #         for e_id in incident_edges:  # For every incident edge e { O(d(v)) }
    #             if (verbose):
    #                 print('incident edge: ', e_id)
    #             if e_id in self.e_indices:
    #                 e = self.get_edge_byindex(e_id)
    #                 for u in e:  # For every node u \in e { O(|e|)}
    #                     if u != v:  # If u != v
    #                         if (verbose):
    #                             print('visiting ', u)
    #                         if not visited_dict.get(u, False):  # { O(1) }
    #                             visited_dict[u] = True  # { O(1) }
    #                             # self.nbr[u].remove(v) # v is no more a neighbour of u, because u will be deleted
    #                             # if len(self.inc_dict[u]) <= 1: # If u has no edge except e incident on it {O (d(nbr(v)))}
    #                             # So a neighbour of the vertex we want to delete has degree = 1
    #                             if self.degree_dict[u] <= 1:
    #                                 if (verbose):
    #                                     print('removing incident edge entry of ', u)
    #                                 del self.inc_dict[u]  # Remove u, because its deg became 0. { O(1) }
    #
    #                             else:
    #                                 if (verbose):
    #                                     print('removing ', e_id, ' from ', u, ' \'s incident edges')
    #                                 self.inc_dict[u].remove(
    #                                     e_id)  # Otherwise, update the incidence entry of u { O (1) }
    #
    #                             self.degree_dict[u] -= 1  # degree of u decreases by 1
    #                             if self.degree_dict[u] == 0:
    #                                 if (verbose):
    #                                     print(u, ' : a 0-degree vertex')
    #                                 # isolated_vertices.add(u) # keep track of isolated vertices.
    #                                 # self.nbr[v].remove(u) # keep track of non-isolated vertices.
    #                                 del self.nbr[u]
    #                                 del self.len_nbr_dict[u]
    #                                 del self.degree_dict[u]
    #                 set_e = set(e)
    #                 for u in self.nbr[v]:
    #                     if u in self.nbr:
    #                         self.nbr[u] -= set_e  # Remove e from nbrs of v.
    #                     if u in self.len_nbr_dict:
    #                         self.len_nbr_dict[u] = len(self.nbr.get(u, set()))
    #
    #                 # for u in isolated_vertices:
    #                 #     del self.degree_dict[u]
    #
    #                 if (verbose):
    #                     print('deleting ', e_id)
    #                 del self.e_indices[e_id]
    #
    #         del self.inc_dict[v]
    #         del self.nbr[v]
    #         del self.len_nbr_dict[v]
    #         del self.degree_dict[v]

    # def strong_subgraph(self, vertex_list):
    #     """ returns: Hypergraph object. """
    #     H = Hypergraph()
    #     e_indices = {}  # (position, edge_size) of edge e in e_nodes list
    #     e_nodes = []  # flattened edge list
    #     inc_dict = {}
    #     H.i = 0
    #
    #     # print('inc_dict: ',self.inc_dict.items())
    #     # print('e_indices: ',self.e_indices.items())
    #     # print('e_nodes: ',self.e_nodes)
    #
    #     for v in vertex_list:
    #         for e_id in self.inc_dict.get(v, []):
    #             if e_id not in e_indices:
    #                 e_start, e_end = self.e_indices[e_id]
    #                 flag = True
    #                 e = self.e_nodes[e_start:e_end]
    #                 for u in e:
    #                     if u not in vertex_list:
    #                         flag = False
    #                         break
    #                 if flag:
    #                     e_nodes += e
    #                     _lene = (e_end - e_start)
    #                     e_indices[e_id] = (H.i, H.i + _lene)
    #                     inc_dict[v] = inc_dict.get(v, []) + [e_id]
    #                     H.i += _lene
    #             else:
    #                 inc_dict[v] = inc_dict.get(v, []) + [e_id]
    #
    #     # print('After: ','inc_dict = ',inc_dict.items(),'\n','e_indices = ',e_indices,'\n',' e_nodes = ',e_nodes)
    #     H.e_nodes = e_nodes
    #     H.inc_dict = inc_dict
    #     H.e_indices = e_indices
    #     return H
    #
    # def strong_subgraph2(self, vertex_list):
    #     """
    #     Another implementation of strong_subgraph
    #     returns: Hypergraph object.
    #     """
    #     H = Hypergraph()
    #     e_indices = {}  # (position, edge_size) of edge e in e_nodes list
    #     e_nodes = []  # flattened edge list
    #     inc_dict = {}
    #     H.i = 0
    #
    #     # print('inc_dict: ',self.inc_dict.items())
    #     # print('e_indices: ',self.e_indices.items())
    #     # print('e_nodes: ',self.e_nodes)
    #
    #     # for v in vertex_list:
    #     #     for e_id in self.inc_dict.get(v,[]):
    #     #         if e_id not in e_indices:
    #     #             e_start, e_end = self.e_indices[e_id]
    #     #             flag = True
    #     #             e = self.e_nodes[e_start:e_end]
    #     #             for u in e:
    #     #                 if u not in vertex_list:
    #     #                     flag = False
    #     #                     break
    #     #             if flag:
    #     #                 e_nodes += e
    #     #                 _lene = (e_end - e_start)
    #     #                 e_indices[e_id] = (H.i, H.i+ _lene)
    #     #                 inc_dict[v] = inc_dict.get(v,[])+[e_id]
    #     #                 H.i+= _lene
    #     #         else:
    #     #             inc_dict[v] = inc_dict.get(v,[])+[e_id]
    #
    #     for e_id in self.e_indices:
    #         flag = True
    #         e = self.get_edge_byindex(e_id)
    #         for u in e:
    #             if u not in vertex_list:
    #                 flag = False
    #                 break
    #         if flag:
    #             e_start, e_end = self.e_indices[e_id]
    #             e_nodes += e
    #             _lene = (e_end - e_start)
    #             e_indices[e_id] = (H.i, H.i + _lene)
    #             H.i += _lene
    #
    #     for e_id in e_indices:
    #         e = self.get_edge_byindex(e_id)
    #         for u in e:
    #             inc_dict[u] = inc_dict.get(u, []) + [e_id]
    #
    #     # print('After: ','inc_dict = ',inc_dict.items(),'\n','e_indices = ',e_indices,'\n',' e_nodes = ',e_nodes)
    #     H.e_nodes = e_nodes
    #     H.inc_dict = inc_dict
    #     H.e_indices = e_indices
    #     return H

        """ returns: Hypergraph object. """
        H = Hypergraph()
        e_indices = {}  # (position, edge_size) of edge e in e_nodes list
        e_nodes = []  # flattened edge list
        inc_dict = {}
        H.i = 0
    
        # print('inc_dict: ',self.inc_dict.items())
        # print('e_indices: ',self.e_indices.items())
        # print('e_nodes: ',self.e_nodes)
    
        for e_id in self.e_indices:
            e = self.get_edge_byindex(e_id)
            flag = True 
            for u in e:
                if u not in vertex_list:
                    flag = False
                    break
            if flag:
                e_nodes += e
                _lene = len(e)
                e_indices[e_id] = (H.i, H.i + _lene)
                for v in e:
                    if v not in inc_dict:
                        inc_dict[v] = set()
                    if e_id not in inc_dict[v]:
                        inc_dict[v].add(e_id)
                        H.degree_dict[v] = H.degree_dict.get(v,0) + 1
                H.i += _lene
            # else:
            #     inc_dict[v] = inc_dict.get(v, []) + [e_id]
    
        # print('After: ','inc_dict = ',inc_dict.items(),'\n','e_indices = ',e_indices,'\n',' e_nodes = ',e_nodes)
        H.e_nodes = e_nodes
        H.inc_dict = inc_dict
        H.e_indices = e_indices
        return H

    def get_hnx_format(self):
        import sys

        sys.path.append("HyperNetX")
        import hypernetx as hnx
        
        _tempH = {}
        for e_id in self.e_indices.keys():
            _tempH[e_id] = self.get_edge_byindex(e_id)
        return hnx.Hypergraph(_tempH)

    def get_clique_graph(self):
        binary_edges = set()
        for e_id in self.e_indices.keys():
            e = self.get_edge_byindex(e_id)
            # print(e)
            if len(e) < 2:
                continue 
            elif len(e) == 2:
                y = tuple(sorted(e))
                if y not in binary_edges:
                    binary_edges.add(y)
                # else:
                #     print('not adding ', y)
            else:
                for x in itertools.combinations(e,2):
                    y = tuple(sorted(x))
                    if y not in binary_edges:
                        binary_edges.add(y)
                    # else:
                    #     print('(comb) not adding ', y)
        scenes = {}
        for i, edge in enumerate(binary_edges):
            scenes[i] = list(edge)
        # print(scenes)
        return Hypergraph(scenes)

    # def weak_subgraph(self, vertex_list):
    #     """ returns: Hypergraph object. """
    #     pass

    def get_N(self):
        """ Return num of vertices """
        return len(self.inc_dict)
    
    def get_M(self):
        """ Return num of edges """
        return len(self.e_indices)
    
    def get_degree_sequence(self):
        """ Return the degree sequence in descending order """
        degs = []
        for v in self.degree_dict:
            degs.append(self.degree_dict[v])
        return sorted(degs,reverse = True)

    def get_degree_distr(self):
        """ Return the degree distribution """
        degs = {}
        N = self.get_N()
        for d in self.degree_dict.values():
            degs[d] = degs.get(d,0)+ (1.0/N)
        return sorted(degs.items(),reverse = True)
        # return sorted(degs.items(), reverse = True, key = lambda x: x[1])

    def get_dim_sequence(self):
        """ Return the dimension sequence in descending order """
        dims = []
        for e_start,e_end in self.e_indices.values():
            dims.append(e_end - e_start)
        return sorted(dims,reverse = True)

    def get_dim_distr(self):
        """ Return the dimension distribution """
        assert isinstance(H, hnx.Hypergraph)
        dims = {}
        M = self.get_M()
        for _dim in self.get_dim_sequence():
            dims[_dim] = dims.get(_dim,0)+ (1.0/M)
        return sorted(dims.items(),reverse = True)

    def get_nbr_sequence(self):
        """ Return the sequence nbrhood sizes  in descending order """
        
        return sorted(self.init_nbrsize.values(),reverse = True)

    def get_nbr_distr(self):
        """ Return the distribution of nbr sizes  """
        nbrs = {}
        N = self.get_N()
        for nbr in self.init_nbrsize.values():
            nbrs[nbr] = nbrs.get(nbr,0) + (1.0/N)
        return sorted(nbrs,reverse = True)

    def get_degree_stats(self):
        """ Return the stats of degrees. """
        import pandas as pd
        deg_seq = self.get_degree_sequence()
        stat = {'mean': None, 'max': None, 'min': None, '25%': None, '50%': None, '75%': None, 'std': None}
        _temp = pd.Series(deg_seq).describe()
        stat['mean'] = _temp['mean']
        stat['std'] = _temp['std']
        stat['min'] = _temp['min']
        stat['max'] = _temp['max']
        stat['25%'] = _temp['25%']
        stat['50%'] = _temp['50%']
        stat['75%'] = _temp['75%']
        return stat

    def get_dim_stats(self):
        """ Return the stats of dimensions. """
        import pandas as pd
        dim_seq = self.get_dim_sequence()
        stat = {'mean': None, 'max': None, 'min': None, '25%': None, '50%': None, '75%': None, 'std': None}
        _temp = pd.Series(dim_seq).describe()
        stat['mean'] = _temp['mean']
        stat['std'] = _temp['std']
        stat['min'] = _temp['min']
        stat['max'] = _temp['max']
        stat['25%'] = _temp['25%']
        stat['50%'] = _temp['50%']
        stat['75%'] = _temp['75%']
        return stat

    def get_nbr_stats(self):
        """ Return the stats of neighbourhoods. """
        import pandas as pd
        nbr_seq = self.get_nbr_sequence()
        stat = {'mean': None, 'max': None, 'min': None, '25%': None, '50%': None, '75%': None, 'std': None}
        _temp = pd.Series(nbr_seq).describe()
        stat['mean'] = _temp['mean']
        stat['std'] = _temp['std']
        stat['min'] = _temp['min']
        stat['max'] = _temp['max']
        stat['25%'] = _temp['25%']
        stat['50%'] = _temp['50%']
        stat['75%'] = _temp['75%']
        return stat

    def WriteDegreeDist(self, filename = "deg_dist"):
        import pickle
        with open(filename+'.pickle','wb') as wf:
            dict_tosave = self.get_degree_distr()
            pickle.dump(dict_tosave, wf, protocol=pickle.HIGHEST_PROTOCOL)

    def __str__(self):
        return ",".join([str(i) for i in self.edge_iterator()])

    def initBiparite(self):
        self.vtoEdgeId = {} # v1 => [e1,e2,...]
        self.eidtov = {} # e1 => [v1,v2,..]
        for eid, edge in self.edge_eid_iterator():
            self.eidtov[eid] = list(edge)
            for v in edge:
                if v not in self.vtoEdgeId:
                    self.vtoEdgeId[v] = [] 
                self.vtoEdgeId[v].append(eid)

    def bipartiteDeleteu(self, u):
        for eid in self.vtoEdgeId[u]:
            if eid in self.eidtov:
                if u in self.eidtov[eid]:
                    self.eidtov[eid].remove(u)
        del self.vtoEdgeId[u]

    def getDist2nbr(self, u):
        s = set()
        for eid in self.vtoEdgeId.get(u,[]):
            for v in self.eidtov.get(eid,[]):
                if v!=u:
                    s.add(v)
        return sorted(list(s)) 