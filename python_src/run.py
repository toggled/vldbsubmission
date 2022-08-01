# import sys
# sys.path.append("HyperNetX")
# import matplotlib.pyplot as plt
import networkx as nx
# import hypernetx as hnx
# from hgDecompose.hgDecompose import HGDecompose
# from hgDecompose.utils import get_hg_hnx
# from hgDecompose.newhgDecompose import HGDecompose
from hgDecompose.optimizedhgDecompose import HGDecompose
from hgDecompose.utils import get_hg,check_connectivity,writeHypergraphHg
from hgDecompose.influence_propagation import propagate_for_all_vertices, propagate_for_random_seeds, run_intervention_exp2,run_intervention_exp2_explain,run_intervention_exp2_explain_splen
import argparse
import pandas as pd
import pickle
import os
from copy import deepcopy

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--thread", help="index of thread", default=-1, type=int)
parser.add_argument("-d", "--dataset", type=str, default="default")
parser.add_argument("-a", "--algo", type=str, default="naive_nbr")
parser.add_argument("-v", "--verbose", action='store_true')
parser.add_argument("-s", "--param_s", help="parameter for improve2_nbr", default=1, type=int)
parser.add_argument("--iterations", help="number of iterations", default=1, type=int)
parser.add_argument("-nt", "--nthreads", help="number of threads for improve3_nbr", default=4, type=int)
parser.add_argument("--sis", action='store_true')
parser.add_argument("--sir", action='store_true')
parser.add_argument("--sir_exp2", action='store_true')
parser.add_argument("--sir_exp3", action='store_true') # intervention
parser.add_argument("--sir_exp3_explanation", action = 'store_true')
parser.add_argument("--sir_exp3_explanation_splen", action = 'store_true')
parser.add_argument("-p", "--prob", help="parameter for Probability", default= 0.3, type=float)
parser.add_argument("-g", "--gamma", help="parameter for Probability", default= 0.01, type=float)

args = parser.parse_args()

# Pandemic propagation
if(args.sir or args.sir_exp2 or args.sir_exp3 or args.sir_exp3_explanation or args.sir_exp3_explanation_splen):

    input_H = get_hg(args.dataset)

    H = deepcopy(input_H)
    assert H is not None

    # Loading/saving to file
    os.system("mkdir -p tests/tmp")
    fname = "tests/tmp/" + args.dataset + "_" + args.algo + ".pkl"
    if(not os.path.isfile(fname)):
        hgDecompose = HGDecompose()
        if(args.algo == "naive_nbr"):
            hgDecompose.naiveNBR(input_H, verbose=args.verbose)
        elif(args.algo == "naive_degree"):
            hgDecompose.naiveDeg(input_H, verbose=args.verbose)
        elif(args.algo == "graph_core"):
            G = H.get_clique_graph()
            nx_G = nx.Graph()
            # print("N: ",G.get_N())
            # print("M: ",G.get_M())
            # hgDecompose.naiveDeg(G, verbose=args.verbose)
            for e in G.edge_iterator():
                nx_G.add_edge(e[0],e[1])
            hgDecompose.core = nx.core_number(nx_G)
        else:

            raise RuntimeError(args.algo + " is not defined or implemented yet")

        core_base = hgDecompose.core
        # print(core_base)

        # dump file
        with open(fname, 'wb') as handle:
            pickle.dump(hgDecompose, handle, protocol= 4)


    else:
        # print("Retrieving saved file")
        with open(fname, 'rb') as handle:
            hgDecompose = pickle.load(handle)
            core_base = hgDecompose.core
    
    # quit()
    # print(core_base)
    entry = {}
    entry['dataset'] = args.dataset
    entry['p'] = float(args.prob)
    entry['algo'] = args.algo
    assert not (args.sir and args.sir_exp2)
    if(args.sir):
        entry['exp'] = "sir"
    elif(args.sir_exp2):
        entry['exp'] = "sir_exp2"
    elif(args.sir_exp3):
        entry['exp'] = "sir_exp3"
    elif(args.sir_exp3_explanation):
        entry['exp'] = 'sir_exp3_explanation'
    elif(args.sir_exp3_explanation_splen):
        entry['exp'] = 'sir_exp3_explanation_splen'
    else:
        raise NotImplementedError()

    entry['result'] = None
    entry['timestep_results'] = None  
    entry['intervention_results'] = None 

    
    
    if(args.sir):
        entry['result'] = propagate_for_all_vertices(H, core_base, p = float(args.prob), verbose=args.verbose)
        result = pd.DataFrame()
        result = result.append(entry, ignore_index=True)
        result.to_csv('../output/propagation_result.csv', header=False,
                            index=False, mode='a')
    elif(args.sir_exp2):
        entry['timestep_results'] = propagate_for_random_seeds(H, core_base, p = float(args.prob), verbose=args.verbose)
        result = pd.DataFrame()
        result = result.append(entry, ignore_index=True)
        result.to_csv('../output/propagation_result.csv', header=False,
                            index=False, mode='a')
    elif(args.sir_exp3):
        # entry['result'], entry['timestep_results'] = propagate_for_random_seeds(H, core_base, p = float(args.prob), verbose=args.verbose)
        # entry['intervention_results'] = run_intervention_exp(H, core_base, p = float(args.prob),verbose = args.verbose)
        # entry['intervention_results'] = run_intervention_exp2(args.dataset+"_"+args.algo, original_n = len(H.nodes()), p = float(args.prob),verbose = args.verbose)
        entry['intervention_results'] = run_intervention_exp2(args.dataset+"_"+args.algo, original_n = None, p = float(args.prob),verbose = args.verbose)
        result = pd.DataFrame()
        result = result.append(entry, ignore_index=True)
        result.to_csv('../output/propagation_result_recursive_delinner_'+args.dataset+"_"+args.algo+'3.csv', header=False,
                         index=False, mode='a')
    elif(args.sir_exp3_explanation):
        run_intervention_exp2_explain(args.dataset+"_"+args.algo,verbose = args.verbose)
        quit()
    elif(args.sir_exp3_explanation_splen):
        run_intervention_exp2_explain_splen(args.dataset+"_"+args.algo,verbose = args.verbose)
        quit()
    result = pd.DataFrame()
    result = result.append(entry, ignore_index=True)
    if(args.verbose): 
        print(entry)
        print("\n")
        print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))

    
    print(result)
    # result.to_csv('data/output/propagation_result_exp3.csv', header=False,
    #                         index=False, mode='a')
    # result.to_csv('data/output/propagation_result_topk_exp3.csv', header=False,
    #                         index=False, mode='a')
    print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))
    # result.to_csv('data/output/propagation_result_topkpercent_exp3.csv', header=False,
    #                     index=False, mode='a')
    # result.to_csv('data/output/propagation_result_recursive_delinner_'+args.dataset+"_"+args.algo+'3.csv', header=False,
    #                      index=False, mode='a')
    quit()
