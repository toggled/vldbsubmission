# import sys
# sys.path.append("HyperNetX")
# import matplotlib.pyplot as plt
import networkx as nx
# import hypernetx as hnx
# from hgDecompose.hgDecompose import HGDecompose
# from hgDecompose.utils import get_hg_hnx
# from hgDecompose.newhgDecompose import HGDecompose
from hgDecompose.optimizedhgDecompose import HGDecompose
from hgDecompose.utils import get_hg, check_connectivity, writeHypergraphHg
# from hgDecompose.influence_propagation import propagate_for_all_vertices, propagate_for_random_seeds, run_intervention_exp2,run_intervention_exp2_explain,run_intervention_exp2_explain_splen, propagate_for_all_vertices_for_kd
from hgDecompose.influence_propagation_new import exp_9a, exp_kd, propagate_for_all_vertices, propagate_for_all_vertices_for_kd, run_intervention_exp2

import argparse
import pandas as pd
import pickle
import os
from copy import deepcopy

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="default")
parser.add_argument("-a", "--algo", type=str, default="naive_nbr")
parser.add_argument("-v", "--verbose", action='store_true')
parser.add_argument(
    "--seed_size", help="number of seed nodes", default=100, type=int)
parser.add_argument("--max_propagation_time",
                    help="number of iterations in sir", default=100, type=int)
parser.add_argument(
    "--iterations", help="number of iterations", default=1, type=int)
parser.add_argument("-sir", "--sir", default=0, type=int)
parser.add_argument("--sir_9a", action='store_true')
parser.add_argument("-sir_kd", "--sir_kd", action='store_true')
parser.add_argument("-sir_dk", "--sir_dk", action='store_true')
parser.add_argument("-sir_exp2", "--sir_exp2", default=0, type=int)
parser.add_argument("-sir_exp3", "--sir_exp3",
                    default=0, type=int)  # intervention
parser.add_argument("-n", "--num_delete", type=int,
                    default=-1, help="how many vertices are deleted")
parser.add_argument("-sir_exp3_explanation",
                    "--sir_exp3_explanation", default=0, type=int)
parser.add_argument("-sir_exp3_explanation_splen",
                    "--sir_exp3_explanation_splen", default=0, type=int)
parser.add_argument(
    "-p", "--prob", help="parameter for Probability", default=0.3, type=float)
parser.add_argument(
    "-g", "--gamma", help="parameter for Probability", default=0.01, type=float)
args = parser.parse_args()

# Pandemic propagation

if(False):
    input_H = get_hg(args.dataset)

    H = deepcopy(input_H)
    assert H is not None

# Loading/saving to file
os.system("mkdir -p tests/tmp")
os.system("mkdir -p ../output/")
fname = "tests/tmp/" + args.dataset + "_" + args.algo + ".pkl"

entry = {}
entry['dataset'] = args.dataset
entry['p'] = float(args.prob)
entry['algo'] = args.algo
assert not (args.sir and args.sir_exp2)
if(args.sir):
    entry['exp'] = "sir"
elif(args.sir_9a):
    entry['exp'] = "sir_9a"
elif(args.sir_kd):
    entry['exp'] = "sir_kd"
elif(args.sir_dk):
    entry['exp'] = "sir_dk"
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
entry['max propagation time'] = None


if(args.sir_9a):

    core_data_filename = "sirdata_naheed_vai/core_" + \
        args.algo + "_" + args.dataset + "_h0_" + \
        str(args.num_delete) + ".csv"

    # get core information
    core_base = {}
    with open(core_data_filename) as file:
        for line in file:
            vs = line.strip().split(",")
            assert int(vs[0]) not in core_base
            core_base[int(vs[0])] = int(vs[1])

    neighbor_data_filename = "sirdata_naheed_vai/" + args.algo + \
        "_" + args.dataset + "_h0_" + str(args.num_delete) + ".csv"

    # get neighborhood information
    neighbor = {}
    with open(neighbor_data_filename) as file:
        for line in file:
            vs = line.strip().split(",")
            assert vs[0] not in neighbor
            neighbor[int(vs[0])] = list(map(int, vs[1:]))
            for u in neighbor[int(vs[0])]:
                assert type(u) == int
            assert int(vs[0]) in core_base

    entry['result'] = exp_9a(neighbor, core_base, num_vertex_per_core=args.seed_size,
                             num_iterations=args.max_propagation_time, p=float(args.prob), verbose=args.verbose)
    entry['max propagation time'] = args.max_propagation_time
    entry['seed size'] = args.seed_size
    entry['num delete'] = args.num_delete

    result = pd.DataFrame()
    result = result.append(entry, ignore_index=True)
    result.to_csv('../output/propagation_result_9a.csv', header=False,
                  index=False, mode='a')
elif(args.sir_kd or args.sir_dk):
    assert args.sir_dk != args.sir_kd
    import pandas as pd

    assert args.num_delete == -1
    neighbor_data_filename = "sirdata_naheed_vai/" + args.algo + \
        "_" + args.dataset + "_h0_" + str(args.num_delete) + ".csv"
    # get neighborhood information
    neighbor = {}
    with open(neighbor_data_filename) as file:
        for line in file:
            vs = line.strip().split(",")
            assert vs[0] not in neighbor
            neighbor[int(vs[0])] = list(map(int, vs[1:]))
            for u in neighbor[int(vs[0])]:
                assert type(u) == int
            # assert int(vs[0]) in core_base

    kd_result = pd.read_csv(
        "sirdata_naheed_vai/core_kdcore_" + args.dataset + ".csv", header=None)
    kd_result.columns = ['vertex', 'k', 'd']

    # kd_result.sort_values(by=['k', 'd'], ascending=False, inplace=True)
    # a dictionary, where key = (k, d) and value is a list of vertex belonging to that (k, d) core
    kd = {}
    for key, item in kd_result.groupby(['k', 'd'], as_index=False):
        # print(item[''])
        assert len(item['vertex'].unique()) == item.shape[0]
        kd[(int(key[0]), int(key[1]))] = list(
            map(int, item['vertex'].unique()))
        # print(kd[(int(key[0]), int(key[1]))])
        for v in kd[(int(key[0]), int(key[1]))]:
            assert type(v) == int
            assert v in neighbor
        # quit()

    # import pprint
    # pprint.pprint(kd)
    # quit()

    algo_name = "kd" if args.sir_kd else "dk"
    entry['algo'] = algo_name
    entry['result'] = exp_kd(algo_name,
                             neighbor, kd, num_vertex_per_core=args.seed_size, num_iterations=args.max_propagation_time, p=float(args.prob), verbose=args.verbose)
    entry['max propagation time'] = args.max_propagation_time
    entry['seed size'] = args.seed_size
    entry['num delete'] = args.num_delete
    result = pd.DataFrame()
    result = result.append(entry, ignore_index=True)
    result.to_csv('../output/propagation_result_9a.csv', header=False,
                  index=False, mode='a')
result = pd.DataFrame()
result = result.append(entry, ignore_index=True)
if(args.verbose):
    print(entry)
    print("\n")
    print(
        ", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))


print(result)
print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))
quit()


# if(not os.path.isfile(fname)):
#     hgDecompose = HGDecompose()
#     if(args.algo == "naive_nbr"):
#         hgDecompose.naiveNBR(input_H, verbose=args.verbose)
#     elif(args.algo == "naive_degree"):
#         hgDecompose.naiveDeg(input_H, verbose=args.verbose)
#     elif(args.algo == "graph_core"):
#         G = H.get_clique_graph()
#         nx_G = nx.Graph()
#         # print("N: ",G.get_N())
#         # print("M: ",G.get_M())
#         # hgDecompose.naiveDeg(G, verbose=args.verbose)
#         for e in G.edge_iterator():
#             nx_G.add_edge(e[0], e[1])
#         hgDecompose.core = nx.core_number(nx_G)
#     else:

#         raise RuntimeError(args.algo + " is not defined or implemented yet")

#     core_base = hgDecompose.core
#     # print(core_base)

#     # dump file
#     with open(fname, 'wb') as handle:
#         print('dump: ', fname)
#         pickle.dump(hgDecompose, handle, protocol=4)

# else:
#     print("Retrieving saved file")
#     try:
#         with open(fname, 'rb') as handle:
#             hgDecompose = pickle.load(handle)
#             core_base = hgDecompose.core
#     except:
#         # for dblp dataset
#         with open(fname, 'rb') as handle:
#             core_base = pickle.load(handle)

#     # revise core_base
#     core_base_new = {}
#     for v in core_base:
#         core_base_new[int(v)] = int(core_base[v])
#     core_base = core_base_new


# result.to_csv('data/output/propagation_result_exp3.csv', header=False,
#                         index=False, mode='a')
# result.to_csv('data/output/propagation_result_topk_exp3.csv', header=False,
#                         index=False, mode='a')


# result.to_csv('data/output/propagation_result_topkpercent_exp3.csv', header=False,
#                     index=False, mode='a')
# result.to_csv('data/output/propagation_result_recursive_delinner_'+args.dataset+"_"+args.algo+'3.csv', header=False,
#                      index=False, mode='a')


# elif(args.sir):
#     entry['result'] = propagate_for_all_vertices(
#         neighbor, core_base, num_iterations=args.max_propagation_time, p=float(args.prob), verbose=args.verbose)
#     entry['max propagation time'] = args.max_propagation_time

#     result = pd.DataFrame()
#     result = result.append(entry, ignore_index=True)
#     result.to_csv('../output/propagation_result.csv', header=False,
#                   index=False, mode='a')
# elif(args.sir_exp2):
#     entry['timestep_results'] = propagate_for_random_seeds(
#         H, core_base, p=float(args.prob), verbose=args.verbose)
#     result = pd.DataFrame()
#     result = result.append(entry, ignore_index=True)
#     result.to_csv('../output/propagation_result.csv', header=False,
#                   index=False, mode='a')
# elif(args.sir_exp3):
#     # entry['result'], entry['timestep_results'] = propagate_for_random_seeds(H, core_base, p = float(args.prob), verbose=args.verbose)
#     # entry['intervention_results'] = run_intervention_exp(H, core_base, p = float(args.prob),verbose = args.verbose)
#     # entry['intervention_results'] = run_intervention_exp2(args.dataset+"_"+args.algo, original_n = len(H.nodes()), p = float(args.prob),verbose = args.verbose)
#     entry['intervention_results'] = run_intervention_exp2(args.dataset+"_"+args.algo + "_" + str(
#         args.num_delete), original_n=None, p=float(args.prob), verbose=args.verbose)
#     # import pprint
#     # pprint.pprint(entry['intervention_results'])
#     entry['num delete'] = args.num_delete
#     result = pd.DataFrame()
#     result = result.append(entry, ignore_index=True)
#     result.to_csv('../output/propagation_result_recursive_delinner_'+args.dataset+"_"+args.algo+'3.csv', header=False,
#                   index=False, mode='a')
# elif(args.sir_exp3_explanation):
#     run_intervention_exp2_explain(
#         args.dataset+"_"+args.algo, verbose=args.verbose)
#     quit()
# elif(args.sir_exp3_explanation_splen):
#     run_intervention_exp2_explain_splen(
#         args.dataset+"_"+args.algo, verbose=args.verbose)
#     quit()
