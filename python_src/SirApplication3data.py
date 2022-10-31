from hgDecompose.utils import get_hg, writeHypergraphHg, writeHypergraph, get_random_hg,save_dictascsv
from hgDecompose.Hypergraph import Hypergraph
import random,os
import time 
import argparse,pickle,copy
from hgDecompose.optimizedhgDecompose import HGDecompose
from matplotlib import pyplot as plt
import networkx as nx 
# random.seed(10)
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="default")
parser.add_argument("-a", "--algo", type=str, default="naive_nbr")
parser.add_argument("-n", "--num_delete", type=int, default=20, help="how many vertices are deleted")
parser.add_argument("-l", "--level", help="how many times innermost core is deleted", default=1, type=int)
args = parser.parse_args()


def del_innercore(H, diction):
    max_core = max(diction.values())
    remainder = {}
    deleted = {}
    for k,v in diction.items():
        if v < max_core:
            remainder[k] = v 
        else:
            deleted[k] = v
    deleted_list = list(deleted.keys())

    print(len(deleted_list))
    
    # delete num_delete nodes only
    num_delete_local = args.num_delete
    if(num_delete_local > len(deleted_list)):
        num_delete_local = len(deleted_list)
    sampled_list = random.sample(deleted_list, len(deleted_list) - num_delete_local)
    for k in sampled_list:
        remainder[k] = deleted[k]
    return remainder

def gen_nested_hypergraph():
    # Delete innermost core 10 times, each time operating on the hypergraph from previous iteration.
    pathstring = "sirdata/"
    output = {} 
    os.system('mkdir -p '+pathstring)
    os.system("mkdir -p tests/tmp")
    
    name = args.dataset
    algoname = args.algo
    level = int(args.level)

    input_H = get_hg(name)
    fname = "tests/tmp/" + name + "_" + algoname + ".pkl"
    if(not os.path.isfile(fname)):
        hgDecompose = HGDecompose()
        if(args.algo == "naive_nbr"):
            hgDecompose.naiveNBR(input_H, verbose=False)
        if(args.algo == "naive_degree"):
            hgDecompose.naiveDeg(input_H, verbose=False)

        if(args.algo == "graph_core"):
            G = input_H.get_clique_graph()
            nx_G = nx.Graph()
            # print("N: ",G.get_N())
            # print("M: ",G.get_M())
            # hgDecompose.naiveDeg(G, verbose=args.verbose)
            for e in G.edge_iterator():
                nx_G.add_edge(e[0],e[1])
            hgDecompose.core = nx.core_number(nx_G)
        core_base = hgDecompose.core
        
    else:
        with open(fname, 'rb') as handle:
            hgDecompose = pickle.load(handle)
            core_base = hgDecompose.core

    output[0] = {}
    output[0]['H'] = input_H
    output[0]['core'] = core_base 
    # output_core_fname = "tests/tmp/" + name + "_" + algoname + "_"+str(0)+".csv"
    # output_hg_fname = "tests/tmp/" + name + "_" + algoname + "_"+str(0)+".hyp"
    # writeHypergraphHg(output[0]['H'],output_hg_fname)
    # save_dictascsv(output[0]['core'],output_core_fname)

    remainder_vertices = del_innercore(input_H, core_base)
    for i in range(1,level):
        print('i: ',i)
        output[i] = {}
        input_H = input_H.strong_subgraph(remainder_vertices)
        _edgedict = {}
        for eid,e in input_H.edge_eid_iterator():
            _edgedict[eid] = e 
        output[i]['H'] = Hypergraph(_edgedict)
        hgDecompose = HGDecompose()
        if(args.algo == "naive_nbr"):
            hgDecompose.naiveNBR(copy.deepcopy(output[i]['H']), verbose=False)
        if(args.algo == "naive_degree"):
            hgDecompose.naiveDeg(copy.deepcopy(output[i]['H']), verbose=False)
        if(args.algo == "graph_core"):
            G = copy.deepcopy(output[i]['H']).get_clique_graph()
            nx_G = nx.Graph()
            # print("N: ",G.get_N())
            # print("M: ",G.get_M())
            # hgDecompose.naiveDeg(G, verbose=args.verbose)
            for e in G.edge_iterator():
                nx_G.add_edge(e[0],e[1])
            hgDecompose.core = nx.core_number(nx_G)
        core_base = hgDecompose.core
        output[i]['core'] = core_base
        remainder_vertices = del_innercore(output[i]['H'], core_base)

        # output_core_fname = "tests/tmp/" + name + "_" + algoname + "_"+str(i)+".csv"
        # output_hg_fname = "tests/tmp/" + name + "_" + algoname + "_"+str(i)+".hyp"
        # writeHypergraphHg(output[i]['H'],output_hg_fname)
        # save_dictascsv(output[i]['core'],output_core_fname)
            
    with open(os.path.join(pathstring,name+'_'+algoname + "_" + str(args.num_delete) +'.pkl'), 'wb') as handle:
            pickle.dump(output, handle, protocol= 4)


gen_nested_hypergraph()

# pathstring = "data/datasets/sirdata/"
# args = parser.parse_args()
# name = args.dataset
# algoname = args.algo
# level = int(args.level)
# with open(os.path.join(pathstring,name+'_'+algoname+'.pkl'), 'rb') as handle:
# # with open(os.path.join(pathstring,name+'_'+algoname+'.pkl'), 'rb') as handle:
#     output = pickle.load(handle)
# print(output.keys())
# print('H0: ',output[0])
# print('H1: ',output[1])
