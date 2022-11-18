# from hgDecompose.optimizedhgDecompose import HGDecompose
from numpy.core.fromnumeric import mean
import pandas as pd
from hgDecompose.Hypergraph import Hypergraph
import random
import heapq
from hgDecompose.heapdict import heapdict
from disjoint_set import DisjointSet
import pickle
import math 

def strong_subgraph(H, vertex_set):
    import sys

    sys.path.append("HyperNetX")
    import hypernetx as hnx
    """
    Returns the strong sub-hypergraph of H induced by vertex_set
    Parameters
    ----------
    H: Hypernetx Hypergraph
    vertex_set: List/set of vertex label

    Returns
    -------
    Hypernetx Hypergraph object
        The strong sub-hypergraph induced by vertex_set
    """
    assert isinstance(H, hnx.Hypergraph)
    if not isinstance(vertex_set, set):
        X = set(vertex_set)
    else:
        X = vertex_set
    _tempdict = {}  # dictionary for induced edges
    for e_id, e_i in H.incidence_dict.items():
        set_e = set(e_i)
        if set_e.issubset(X):  # If an edge of the original hypergraph is a subset of the vertex set, add it
            _tempdict[e_id] = e_i
    return hnx.Hypergraph(_tempdict)


def get_number_of_nbrs(H, u):
    """
        Returns the number of neighbours of u in hypergraph H
        Parameters
        ----------
        H: Hypernetx Hypergraph
        u: a vertex label

        Returns
        -------
        Integer
            The number of neighbours of u in H
    """
    nbrs = H.neighbors(u)
    if nbrs is None:  # u is not in H
        return 0
    return len(nbrs)


def get_degree(H, u):
    degree = 0
    try:
        degree = H.degree(u)
    except Exception as e:
        # print(e)
        pass
    
    return degree

def get_nbrs(H, u):
    import sys

    sys.path.append("HyperNetX")
    import hypernetx as hnx
    """
        Returns the neighbours of u in hypergraph H
        Parameters
        ----------
        H: Hypernetx Hypergraph
        u: a vertex label

        Returns
        -------
        List
            The neighbours of u in H. [] if u is not in H.
    """
    nbrs = H.neighbors(u)
    if nbrs is None:  # u is not in H
        return []
    return nbrs

def get_hg(dataset):
    H = None
    if(dataset == "default"):
        dic = {
            0: ('FN', 'TH'),
            1: ('TH', 'JV'),
            2: ('BM', 'FN', 'JA'),
            3: ('JV', 'JU', 'CH', 'BM'),
            4: ('JU', 'CH', 'BR', 'CN', 'CC', 'JV', 'BM'),
            5: ('TH', 'GP'),
            6: ('GP', 'MP'),
            7: ('MA', 'GP')
        }

        H = Hypergraph(dic)

    elif(dataset in ['enron', "syn", "klay"]):
        # file location
        dataset_to_filename = {
            # real
            "enron" : "../data/datasets/real/Enron.hyp",
            "klay" : "../data/klay_s.hyp"
        }
        

        # # split by
        # dataset_to_split = {
        #     "enron" : ",",
        #     "congress" : ",",
        #     "contact" : ",",
        #     "dblp": ",",
        #     "amazon": ",",
        #     "syn" : ",",
        #     "bin_1" : ",",
        #     "bin_2" : ",",
        #     "bin_4" : ",",
        #     "bin_5" : ",",

        # }

        
        dic = {}
        # read from file
        with open(dataset_to_filename[dataset]) as f:
            idx = 0
            for line in f:
                edge = tuple(line[:-1].split(','))
                dic[idx] = edge
                idx+=1
                # if idx%10000 == 0:
                #     print(idx)

        H = Hypergraph(dic)

    else:
        raise RuntimeError(dataset + " is not defined or implemented yet")


    return H

def get_random_hg(n = 10, m = 5, edge_size_ub = None, seed = 1):
    """ 
    Returns a random hypergraph with n vertices and m edges. 
    Generate V = {1,2,..,n} 
    Generate E = Take m randomly chosen subsets of V.
    If edge_size_ub is not None, assume every hyperedge can have at most edge_size_ub vertices in it.
    """
    random.seed(seed)
    V = set(range(1, n+1))
    Edict = {}
    if m == 1:
        if n <= edge_size_ub:
            e = tuple([str(i) for i in range(1,n+1)])
            Edict[m] = e
            return Hypergraph(Edict)
    while m:
        if edge_size_ub is None:
            edge_sz = random.randint(2,n) # 2 because we do not want singletons
        else:
            edge_sz = random.randint(2, edge_size_ub) # 2 because we do not want singletons
        
        e = random.sample(V, edge_sz)
        Edict[m] = tuple([str(v) for v in sorted(list(e))])
        m = m-1
    return Hypergraph(Edict)

def get_random_graph(n = 10, m = 5, seed = 1):
    """ 
    Returns a random hypergraph with n vertices and m edges. 
    Generate V = {1,2,..,n} 
    Generate E = Take m randomly chosen subsets of V.
    If edge_size_ub is not None, assume every hyperedge can have at most edge_size_ub vertices in it.
    """
    random.seed(seed)
    V = set(range(1, n+1))
    Edict = {}
    exists = {}
    while m:
        e = random.sample(V, 2)
        if e[1]<e[0]:
            e[0],e[1] = e[1],e[0]
        if not exists.get((e[0],e[1]),False):
            exists[(e[0],e[1])] = True
            Edict[m] = (e[0],e[1])
            m = m-1
    return Hypergraph(Edict)

def get_basic_hg(dataset):
    from hgDecompose.BasicHypergraph import Hypergraph as HypergraphBasic
    H = None
    if(dataset == "default"):
        dic = {
            0: ('FN', 'TH'),
            1: ('TH', 'JV'),
            2: ('BM', 'FN', 'JA'),
            3: ('JV', 'JU', 'CH', 'BM'),
            4: ('JU', 'CH', 'BR', 'CN', 'CC', 'JV', 'BM'),
            5: ('TH', 'GP'),
            6: ('GP', 'MP'),
            7: ('MA', 'GP')
        }

        H = HypergraphBasic(dic)

    elif(dataset in ['enron', "syn", "bin_1", "bin_2", "bin_4", "bin_5", "4_sim", "5_sim", "pref", "pref_20000","pref_40000","pref_60000","pref_80000","pref_100000", "congress", "contact","dblp", "amazon"]):

        # file location
        dataset_to_filename = {
            # real
            "enron" : "data/datasets/real/Enron.hyp",
            "congress" : "data/datasets/real/congress-bills.hyp",
            "contact" : "data/datasets/real/contact-primary-school.hyp",
            "dblp": "data/datasets/real/DBLP.hyp",
            "amazon": "data/datasets/real/amazon-reviews.hyp",

            # synthetic
            "syn" : "data/datasets/synthetic/syn.hyp",
            "bin_1" : "data/datasets/synthetic/binomial_5_100_4_0.200000_sample_1_iter_1.txt",
            "bin_2" : "data/datasets/synthetic/binomial_5_500_4_0.200000_sample_2_iter_1.txt",
            "bin_4" : "data/datasets/synthetic/binomial_5_100_3_0.200000_sample_4_iter_1.txt",
            "bin_5" : "data/datasets/synthetic/binomial_5_500_3_0.200000_sample_5_iter_1.txt",
            "4_sim": "data/datasets/synthetic/4simplex.hyp",
            "5_sim": "data/datasets/synthetic/5simplex.hyp",
            "pref": "data/datasets/synthetic/pref_1000000_3_1.hyp",
            "pref_20000": "data/datasets/synthetic/pref_20000_3_1_simple.hyp",
             "pref_40000": "data/datasets/synthetic/pref_40000_3_1_simple.hyp",
             "pref_60000": "data/datasets/synthetic/pref_60000_3_1_simple.hyp",
             "pref_80000": "data/datasets/synthetic/pref_80000_3_1_simple.hyp",
             "pref_100000": "data/datasets/synthetic/pref_100000_3_1_simple.hyp"
        }

        
        # # split by
        # dataset_to_split = {
        #     "enron" : ",",
        #     "congress" : ",",
        #     "contact" : ",",
        #     "dblp": ",",
        #     "amazon": ",",
        #     "syn" : ",",
        #     "bin_1" : ",",
        #     "bin_2" : ",",
        #     "bin_4" : ",",
        #     "bin_5" : ",",

        # }

        
        dic = {}
        # read from file
        with open(dataset_to_filename[dataset]) as f:
            idx = 0
            for line in f:
                edge = tuple(line[:-1].split(','))
                dic[idx] = edge
                idx+=1
                # if idx%10000 == 0:
                #     print(idx)

        H = HypergraphBasic(dic)

    else:
        raise RuntimeError(dataset + " is not defined or implemented yet")


    return H

def writeHypergraph(edge_dict, out_file):
    with open(out_file,'w') as wf:
        for edge in edge_dict.values():
            edge_str = ",".join([str(node) for node in edge])
            wf.write(edge_str+"\n")

def writeHypergraphHg(hg, out_file):
    assert isinstance(hg, Hypergraph)
    with open(out_file,'w') as wf:
        for id in sorted(list(hg.e_indices.keys())):
            edge = hg.get_edge_byindex(id)
            edge_str = ",".join([str(node) for node in edge])
            wf.write(edge_str+"\n")

def get_N(H):
    """ Return num of vertices """
    return len(H.nodes)

def get_M(H):
    """ Return num of edges """
    return len(H.edges)

def get_degree_sequence(H):
    """ Return the degree sequence in descending order """
    assert isinstance(H, hnx.Hypergraph)
    degs = []
    for v in H.nodes:
        degs.append(H.degree(v))
    return sorted(degs,reverse = True)

def get_degree_distr(H):
    """ Return the degree distribution """
    assert isinstance(H, hnx.Hypergraph)
    degs = {}
    N = get_N(H)
    for v in H.nodes:
        d = H.degree(v)
        degs[d] = degs.get(d,0)+ (1.0/N)
    return sorted(degs.items(),reverse = True)

def get_dim_sequence(H):
    """ Return the dimension sequence in descending order """
    assert isinstance(H, hnx.Hypergraph)
    dims = []
    for e in H.edges:
        dims.append(H.dim(e)+1)
    return sorted(dims,reverse = True)

def get_dim_distr(H):
    """ Return the dimension distribution """
    assert isinstance(H, hnx.Hypergraph)
    dims = {}
    M = get_M(H)
    for _dim in get_dim_sequence(H):
        dims[_dim] = dims.get(_dim,0)+ (1.0/M)
    return sorted(dims.items(),reverse = True)

def get_nbr_sequence(H):
    """ Return the sequence nbrhood sizes  in descending order """
    assert isinstance(H, hnx.Hypergraph)
    nbrs = []
    for v in H.nodes:
        nbrs.append(get_number_of_nbrs(H,v))
    return sorted(nbrs,reverse = True)

def get_nbr_distr(H):
    """ Return the distribution of nbr sizes  """
    assert isinstance(H, hnx.Hypergraph)
    nbrs = {}
    N = get_N(H)
    for nbr in get_nbr_sequence(H):
        nbrs[nbr] = nbrs.get(nbr,0) + (1.0/N)
    return sorted(nbrs,reverse = True)

def get_degree_stats(H):
    """ Return the stats of degrees. """
    assert isinstance(H, hnx.Hypergraph)
    deg_seq = get_degree_sequence(H)
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

def get_dim_stats(H):
    """ Return the stats of dimensions. """
    assert isinstance(H, hnx.Hypergraph)
    dim_seq = get_dim_sequence(H)
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

def get_nbr_stats(H):
    """ Return the stats of neighbourhoods. """
    assert isinstance(H, hnx.Hypergraph)
    nbr_seq = get_nbr_sequence(H)
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


import numpy as np


def quickselect_median(l, pivot_fn=random.choice):
    if len(l) % 2 == 1:
        return quickselect(l, len(l) // 2, pivot_fn)
    else:
        return int(0.5 * (quickselect(l, len(l) / 2 - 1, pivot_fn) +
                    quickselect(l, len(l) / 2, pivot_fn)))

def quickselect(l, k, pivot_fn):
    """
    Select the kth element in l (0 based)
    :param l: List of numerics
    :param k: Index
    :param pivot_fn: Function to choose a pivot, defaults to random.choice
    :return: The kth element of l
    """
    if len(l) == 1:
        assert k == 0
        return l[0]

    pivot = pivot_fn(l)

    lows = [el for el in l if el < pivot]
    highs = [el for el in l if el > pivot]
    pivots = [el for el in l if el == pivot]

    if k < len(lows):
        return quickselect(lows, k, pivot_fn)
    elif k < len(lows) + len(pivots):
        # We got lucky and guessed the median
        return pivots[0]
    else:
        return quickselect(highs, k - len(lows) - len(pivots), pivot_fn)

# https://towardsdatascience.com/fastest-way-to-calculate-h-index-of-publications-6fd52e381fee
# Expert algorithm derived from wiki https://en.wikipedia.org/wiki/H-index
def operator_H(citations):
    
    citations = np.array(citations)
    n         = citations.shape[0]
    array     = np.arange(1, n+1)
    
    # reverse sorting
    citations = np.sort(citations)[::-1]
    # print(citations)

    # intersection of citations and k
    h_idx = np.max(np.minimum(citations, array)) # inside np.minimum is element-wise
    # print(np.minimum(citations, array))
    return h_idx

def par_operator_H(args):
    node, citations = args
    citations = np.array(citations)
    n         = citations.shape[0]
    array     = np.arange(1, n+1)
    
    # reverse sorting
    citations = np.sort(citations)[::-1]
    # print(citations)

    # intersection of citations and k
    h_idx = np.max(np.minimum(citations, array)) # inside np.minimum is element-wise
    # print(np.minimum(citations, array))
    return (node, h_idx)

def operator_H_new(citations):
    len_citations = len(citations)
    s = [0 for _ in range(len_citations + 1)]
    for i in range(len_citations):
        s[min(len_citations, citations[i])] += 1
    # the i-th index in s  is the count of i's (or the smallest between i and len_citations) in citations

    # print(s)
    sum = 0
    for i in range(len_citations - 1, -1, -1):
        sum += s[i]
        if(sum >= i):
            return i
        
    return 0


def operator_H_quicksort(citations):
        
    len_citations = len(citations)
    median = quickselect_median(citations)
    print(citations, median)
    if(len_citations % 2 == 1):
        if(median == (len_citations - 1) / 2):
            return median
        elif(median > (len_citations - 1) / 2):
            return operator_H_new(citations[:median])
        elif(median < (len_citations - 1) / 2):
            return operator_H_new(citations[median:])
    else:
        # When len is even, we may not identify median correctly. median([1, 2]) = 1 but median([1,1,2,2]) = 1 raises error 
        if(median == len_citations / 2):
            return median
        elif(median > len_citations / 2):
            return operator_H_new(citations[:median])
        elif(median < len_citations / 2):
            return operator_H_new(citations[median:])


def operator_H_priority_queue(citations):

    len_citations = len(citations)
    citations = [-1 * citation for citation in citations] # value is negated to implement max-heap
    heapq.heapify(citations) # O(n)

    h = 0
    for i in range(1, len_citations + 1):
        if(-1 * heapq.heappop(citations) >= i):
            h = i
        else:
            break
    return h



# memory profile
def memory_usage_psutil():
    # return the memory usage in MB
    import psutil, os
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)
    return mem

def check_connectivity(hg):
    ds = DisjointSet()
    for v in hg.nodes():
        if v in ds:
            continue 
        else:
            for e in hg.inc_dict[v]:
                # print(hg.get_edge_byindex(e))
                for u in hg.get_edge_byindex(e):
                    ds.union(v,u)
            # print('-----')
            # print(list(ds.itersets()),'\n------')
    print('Is connected: ', len(list(ds.itersets()))==1)
    # print(list(ds.itersets()))
            
def component_sz(v,hg):
    """ Returns the size of the component v is part of in Hg"""
    ds = DisjointSet()
    queue = [v]
    traversed = {v: True}
    while len(queue):
        v = queue.pop(0)
        for e in hg.inc_dict[v]:
            for u in hg.get_edge_byindex(e):
                ds.union(v,u)
                if u not in traversed:
                    queue.append(u)
                    traversed[u] = True
                if not traversed[u]:
                    queue.append(u)
                    traversed[u] = True
    # print('traversed= ',len(traversed), '|V| = ', len(hg.inc_dict))
    # print(len(list(ds.itersets())))
    # assert len(traversed) == len(hg.inc_dict)
    return len(traversed)

def avg_shortest_pathlen(source, hg, number_of_targets, V, constant_M, verbose = True):
    """ Dijkstra's algorithm on hypergraph """
    # V = list(hg.inc_dict.keys())
    dist = {}
    traversed = {}
    for u in V:
        dist[u] = math.inf
        traversed[u] = False 
    # for u in random.choices(V,k = number_of_targets):
    #     # Compute number of edges in between u~v 
    #     dist[u] = math.inf

    dist[source] = 0
    heap = heapdict()
    heap[source] = 0
    traversed[source] = True 
    while len(heap):
        v, dist_v = heap.popitem()
        for e in hg.inc_dict[v]:
            edge = hg.get_edge_byindex(e)
            for u in edge:
                if not traversed[u]:
                    dist_root =  dist_v + 1    
                    if dist_root < dist[u]:
                        dist[u] = dist_root
                        heap[u] = dist_root
    # if (verbose):
    #     print('shortest paths: ',source,' = ')
    #     print(dist)

    # dists = [ dist[u] for u in random.choices(V,k = number_of_targets) ]
    dists = [0]*number_of_targets
    i = 0

    while i < number_of_targets:
        u = random.choice(V)
        if dist[u] is not math.inf:  
            dists[i] = dist[u]
        else:
            dists[i] = constant_M
        i+=1
    return np.mean(dists)


    # while i<number_of_targets:
    #     u = random.choice(V)
    #     if dist[u] is not math.inf:
    #         dists[i] = dist[u]
    #         i+=1
    # return np.mean(dists)

def save_dict(dict,fname = 'tests/tmp/temp.pkl'):
    """ Save a dictionary """
    print('Saving dictionary to: ',fname)
    with open(fname, 'wb') as handle:
        pickle.dump(dict, handle, protocol= 4) 

def load_dict(fname = 'tests/tmp/temp.pkl'):
    """ Load a dictionary """
    print('Loading dictionary from: ',fname)
    with open(fname, 'rb') as handle:
        dict = pickle.load(handle)
        return dict 

def save_dictascsv(dict,fname = 'tests/tmp/temp.csv'):
    """ Save a dictionary """
    print('Saving dictionary to: ',fname)
    with open(fname, 'w') as handle:
        for k,val in sorted(list(dict.items())):
            handle.write(str(k)+','+str(val)+'\n')