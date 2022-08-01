# import sys
# sys.path.append("HyperNetX")
# import matplotlib.pyplot as plt
import networkx as nx
# import hypernetx as hnx
# from hgDecompose.hgDecompose import HGDecompose
# from hgDecompose.utils import get_hg_hnx
# from hgDecompose.newhgDecompose import HGDecompose
from hgDecompose.optimizedhgDecompose import HGDecompose
from hgDecompose.utils import get_hg, memory_usage_psutil,get_localhg,check_connectivity,writeHypergraphHg
from hgDecompose.influence_propagation import propagate_for_all_vertices, propagate_for_random_seeds, run_intervention_exp2,run_intervention_exp2_explain,run_intervention_exp2_explain_splen
from hgDecompose.sis_propagation import propagateSIS_for_all_vertices
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
parser.add_argument("--con", help="Is connected hypergraph", action='store_true')
parser.add_argument("-p", "--prob", help="parameter for Probability", default= 0.3, type=float)
parser.add_argument("-g", "--gamma", help="parameter for Probability", default= 0.01, type=float)
parser.add_argument('--degD', action='store_true')
parser.add_argument('--volD', action='store_true')
args = parser.parse_args()

if (args.degD):
    input_H = get_hg(args.dataset)

    H = deepcopy(input_H)
    assert H is not None
    os.system("mkdir -p tests/density")
    fname = "tests/density/" + args.dataset + "_deg.hyp"
    hgDecompose = HGDecompose()
    subH,maximal_density = hgDecompose.degree_densest_subh(H,verbose = args.verbose)
    writeHypergraphHg(subH,fname)
    quit() 
if (args.con):
    input_H = get_hg(args.dataset)
    check_connectivity(input_H)
    quit()

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
    elif(args.sir_exp2):
        entry['timestep_results'] = propagate_for_random_seeds(H, core_base, p = float(args.prob), verbose=args.verbose)
    elif(args.sir_exp3):
        # entry['result'], entry['timestep_results'] = propagate_for_random_seeds(H, core_base, p = float(args.prob), verbose=args.verbose)
        # entry['intervention_results'] = run_intervention_exp(H, core_base, p = float(args.prob),verbose = args.verbose)
        # entry['intervention_results'] = run_intervention_exp2(args.dataset+"_"+args.algo, original_n = len(H.nodes()), p = float(args.prob),verbose = args.verbose)
        entry['intervention_results'] = run_intervention_exp2(args.dataset+"_"+args.algo, original_n = None, p = float(args.prob),verbose = args.verbose)
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

    os.system("mkdir -p data/output")
    result.to_csv('data/output/propagation_result.csv', header=False,
                            index=False, mode='a')
    # result.to_csv('data/output/propagation_result_exp3.csv', header=False,
    #                         index=False, mode='a')
    # result.to_csv('data/output/propagation_result_topk_exp3.csv', header=False,
    #                         index=False, mode='a')
    print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))
    # result.to_csv('data/output/propagation_result_topkpercent_exp3.csv', header=False,
    #                     index=False, mode='a')
    # result.to_csv('data/output/propagation_result_recursive_delinner_'+args.dataset+"_"+args.algo+'3.csv', header=False,
    #                     index=False, mode='a')
    quit()

# # Pandemic propagation
# if(args.sis):

#     input_H = get_hg(args.dataset)

#     H = deepcopy(input_H)
#     assert H is not None

#     # Loading/saving to file
#     os.system("mkdir -p tests/tmp")
#     fname = "tests/tmp/" + args.dataset + "_" + args.algo + ".pkl"
#     if(not os.path.isfile(fname)):
#         hgDecompose = HGDecompose()
#         if(args.algo == "naive_nbr"):
#             hgDecompose.naiveNBR(input_H, verbose=args.verbose)
#         elif(args.algo == "naive_degree"):
#             hgDecompose.naiveDeg(input_H, verbose=args.verbose)
#         elif(args.algo == "graph_core"):
#             G = H.get_clique_graph()
#             nx_G = nx.Graph()
#             # print("N: ",G.get_N())
#             # print("M: ",G.get_M())
#             # hgDecompose.naiveDeg(G, verbose=args.verbose)
#             for e in G.edge_iterator():
#                 nx_G.add_edge(e[0],e[1])
#             hgDecompose.core = nx.core_number(nx_G)
#         else:

#             raise RuntimeError(args.algo + " is not defined or implemented yet")

#         core_base = hgDecompose.core
#         # print(core_base)

#         # dump file
#         with open(fname, 'wb') as handle:
#             pickle.dump(hgDecompose, handle, protocol= 4)


#     else:
#         # print("Retrieving saved file")
#         with open(fname, 'rb') as handle:
#             hgDecompose = pickle.load(handle)
#             core_base = hgDecompose.core
    
#     # print(core_base)
#     entry = {}
#     entry['dataset'] = args.dataset
#     entry['p'] = float(args.prob)
#     entry["gamma"] = float(args.gamma)
#     entry['algo'] = args.algo
#     entry['result'] = propagateSIS_for_all_vertices(H, core_base, p = float(args.prob), gamma = float(args.gamma), verbose=args.verbose)

#     result = pd.DataFrame()
#     result = result.append(entry, ignore_index=True)
#     if(args.verbose): 
#         print(entry)
#         print("\n")
#         print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))

#     os.system("mkdir -p data/output")
#     result.to_csv('data/output/sis_propagation_result.csv', header=False,
#                             index=False, mode='a')

#     quit()

# hyper-graph construction
# H = get_hg_hnx(args.dataset)
if args.algo.startswith('opt_local_core'):
    input_H = get_localhg(args.dataset)
else:
    input_H = get_hg(args.dataset)
print("HG construction done!")
assert input_H is not None



for iteration in range(args.iterations):
    H = deepcopy(input_H)
    entry = {}
    entry['algo'] = args.algo
    entry['dataset'] = args.dataset
    entry['num_threads'] = args.nthreads
    # run algo
    hgDecompose = HGDecompose()
    if(args.algo == "naive_nbr"):
        hgDecompose.naiveNBR(H, verbose=args.verbose)

    elif(args.algo == "improved_nbr"):
        hgDecompose.improvedNBR(H, verbose=args.verbose)

    elif(args.algo == "improved_nbr_simple"):
        hgDecompose.improvedNBR_simplified(H, verbose=False)
        
    elif(args.algo == "naive_degree"):
        hgDecompose.naiveDeg(H, verbose=args.verbose)

    elif(args.algo == "graph_core"):
        G = H.get_clique_graph()
        nx_G = nx.Graph()
        # print("N: ",G.get_N())
        # print("M: ",G.get_M())
        # hgDecompose.naiveDeg(G, verbose=args.verbose)
        for e in G.edge_iterator():
            nx_G.add_edge(e[0],e[1])
        hgDecompose.core = nx.core_number(nx_G)

        # hgDecompose.naiveNBR(G, verbose=args.verbose)
        # # for core_num in sorted(hgDecompose.core.values()):
        # #     verify_kcore.verify_subgraph(G, core_num, hgDecompose.core)
        # # print(hgDecompose.core)
        # # G = H.get_clique_graph()
        # # hgDecompose.naiveDeg(G, verbose=args.verbose)
        # print(hgDecompose.core)

    elif(args.algo == "improved2_nbr"):
        assert args.param_s > 0 # Is this assertion valid?
        hgDecompose.improved2NBR(H, s=args.param_s, verbose=args.verbose)

    elif (args.algo == 'par_improved2_nbr'):
        assert args.param_s > 0 # Is this assertion valid?
        hgDecompose.parallel_improved2NBR(H, s=args.param_s, num_threads = args.nthreads, verbose=args.verbose)
    
    elif (args.algo == 'par_improved3_nbr'):
        assert args.param_s > 0 # Is this assertion valid?
        hgDecompose.parallel_improved3NBR(H, s=args.param_s, num_threads = args.nthreads, verbose=args.verbose)
    
    elif(args.algo == "recursive_local_core"):
        hgDecompose.local_core(H, verbose=args.verbose)

    elif(args.algo == "iterative_local_core"):
        hgDecompose.iterative_local_core(H, verbose=args.verbose)

    elif(args.algo == "bst_local_core"):
        hgDecompose.bst_local_core(H, verbose=args.verbose)

    elif(args.algo == "improved_local_core"):
        hgDecompose.improved_local_core(H, verbose=args.verbose)
    
    elif(args.algo == "opt_local_core"):
        # Run local_core algorithm while storing core-correction ammount and other auxiliary information.
        hgDecompose.opt_local_core(H, verbose=args.verbose, store_core_information=True, filename="data/output/"+args.dataset+"_local_core.csv", info_dic={'algo' : args.algo, 'dataset' : args.dataset, 'num_threads' : args.nthreads, 'outer iteration' : iteration})
        
        # Run local_core algorithm without storing other auxiliary information.
        # hgDecompose.opt_local_core(H, verbose=args.verbose, store_core_information=False)

    elif(args.algo == "opt_local_core_fast"):
        # Run local_core algorithm without storing any information except execution time
        hgDecompose.opt_local_core_bare_min(H)

    elif(args.algo == "opt_local_coreI"):
        hgDecompose.opt_local_core_I(H, verbose=args.verbose,store_core_information=True, filename="data/output/"+args.dataset+"_local_coreI.csv", info_dic={'algo' : args.algo, 'dataset' : args.dataset, 'num_threads' : args.nthreads, 'outer iteration' : iteration})

    elif(args.algo == "opt_local_coreII"):
        hgDecompose.opt_local_coreII(H, verbose=args.verbose,store_core_information=True, filename="data/output/"+args.dataset+"_local_coreII.csv", info_dic={'algo' : args.algo, 'dataset' : args.dataset, 'num_threads' : args.nthreads, 'outer iteration' : iteration})

    elif(args.algo == "opt_local_coreIII"):
        hgDecompose.opt_local_coreIII(H, verbose=args.verbose,store_core_information=True, filename="data/output/"+args.dataset+"_local_coreIII.csv", info_dic={'algo' : args.algo, 'dataset' : args.dataset, 'num_threads' : args.nthreads, 'outer iteration' : iteration})

    elif(args.algo == "par_local_core"):
        hgDecompose.par_local_core(H, verbose=args.verbose)

    else:
        raise RuntimeError(args.algo + " is not defined or implemented yet")


    entry['core'] = hgDecompose.core
    entry['param_s'] = args.param_s
    entry['execution time'] = hgDecompose.execution_time
    entry['bucket update time'] = hgDecompose.bucket_update_time
    entry['neighborhood call time'] = hgDecompose.neighborhood_call_time
    entry['degree call time'] = hgDecompose.degree_call_time
    entry['num bucket update'] = hgDecompose.num_bucket_update
    entry['num neighborhood computation'] = hgDecompose.num_neighborhood_computation
    entry['num degree computation'] = hgDecompose.num_degree_computation
    entry['subgraph computation time'] = hgDecompose.subgraph_time
    entry['num subgraph call'] = hgDecompose.num_subgraph_call
    entry['init time'] = hgDecompose.init_time
    entry['outerloop time'] = hgDecompose.loop_time
    entry['total iteration'] = hgDecompose.total_iteration
    entry['inner iteration'] = hgDecompose.inner_iteration
    entry['core_correction time'] = hgDecompose.core_correct_time
    entry['h_index_time'] = hgDecompose.h_index_time
    entry['tau'] = hgDecompose.max_n  # For #iterations vs dataset barplot
    entry['core_correction_volume'] = hgDecompose.core_correctionvol_n #  core_corrections volume per iteration => Ammount of core_correction done. => Relation with runtime
    entry['sum_core_correction_volume'] = hgDecompose.core_correction_volume  # For core_correction volume vs dataset plot
    entry['reduction_in_hhat']  = hgDecompose.reduction_hhat_n  # [ hhat^{n-1} - hhat^{n}, for n \in [1, tau] ] => Convergence plot.

    if(True):
        entry['memory taken'] = memory_usage_psutil()
    # print(entry)
    result = pd.DataFrame()
    result = result.append(entry, ignore_index=True)
    # print('result: ',result['num subgraph call'].values[0],',',result['subgraph computation time'].values[0])
    # print('tolist(): ',result.columns.tolist())
    # print("iter: ",iteration)
    if(args.verbose and iteration==0): 
        # print(entry)
        print("\n")
        print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))

    os.system("mkdir -p data/output")
    # print(result[['dataset','execution time','algo']])
    # continue
    # result.to_csv('data/output/scal_result.csv', header=False,
    #                         index=False, mode='a')
    result.to_csv('data/output/result_imp_lc2.csv', header=False,
                            index=False, mode='a')
    # result.to_csv('data/output/result_protein.csv', header=False,
    #                         index=False, mode='a')
    # result.to_csv('data/output/result_gowalla.csv', header=False,
    #                         index=False, mode='a')
    print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))


    # print(memory_usage_psutil())
    # print(memory_usage_psutil())