import argparse
import pandas as pd
import pickle
import os
from hgDecompose.influence_propagation_new import bfs_bounded
from copy import deepcopy
import random
from tqdm import tqdm
import multiprocessing

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="enron")
parser.add_argument("-a", "--algo", type=str, default="naive_nbr")
parser.add_argument("-v", "--verbose", action='store_true')
parser.add_argument(
    "--seed_size", help="number of seed nodes", default=100, type=int)
parser.add_argument("--max_propagation_time",
                    help="number of iterations in sir", default=100, type=int)
parser.add_argument("--sir_prev_9a", action='store_true')
parser.add_argument("-sir_prev_kd", "--sir_prev_kd", action='store_true')
parser.add_argument("-sir_prev_dk", "--sir_prev_dk", action='store_true')
parser.add_argument(
    "-p", "--prob", help="parameter for Probability", default=0.3, type=float)
parser.add_argument("-n", "--num_delete", type=int,
                    default=-1, help="how many vertices are deleted")
args = parser.parse_args()
os.system("mkdir -p tests/tmp")
os.system("mkdir -p ../output/")


num_rounds = 1


def worker(arg_tuple):

    distinct_core_numbers, core_to_vertex_map, neighbor = arg_tuple

    # propogate for each core
    result = {}
    for c in tqdm(distinct_core_numbers):
        # propagate
        result_single_run = []
        for v in tqdm(core_to_vertex_map[c]):
            result_single_run.append(bfs_bounded(
                neighbor, starting_vertex=v, p=float(args.prob), num_iterations=args.max_propagation_time,  verbose=False))
        result[c] = result_single_run
    return result


if __name__ == '__main__':
    entry = {}
    entry['dataset'] = args.dataset
    entry['p'] = float(args.prob)
    entry['algo'] = args.algo
    if(args.sir_prev_9a):
        entry['exp'] = "sir_prev_9a"
    elif(args.sir_prev_kd):
        entry['exp'] = "sir_prev_kd"
        entry['algo'] = "kd"
    elif(args.sir_prev_dk):
        entry['exp'] = "sir_prev_dk"
        entry['algo'] = "dk"
    else:
        raise NotImplementedError()

    entry['result'] = None
    entry['timestep_results'] = None
    entry['intervention_results'] = None
    entry['max propagation time'] = None


    core_data_filename = "sirdata_naheed_vai/core_" + \
        args.algo + "_" + args.dataset + "_h0_" + \
        str(args.num_delete) + ".csv"

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

    # get core information
    if(args.sir_prev_9a):
        core = {}
        with open(core_data_filename) as file:
            for line in file:
                vs = line.strip().split(",")
                assert int(vs[0]) not in core
                core[int(vs[0])] = int(vs[1])

        core_to_vertex_map = {}
        distinct_core_numbers = []
        for v in core:
            if(core[v] not in core_to_vertex_map):
                core_to_vertex_map[core[v]] = [v]
                distinct_core_numbers.append(core[v])
            else:
                core_to_vertex_map[core[v]].append(v)

    elif(args.sir_prev_kd or args.sir_prev_dk):
        kd_result = pd.read_csv(
            "sirdata_naheed_vai/core_kdcore_" + args.dataset + ".csv", header=None)
        kd_result.columns = ['vertex', 'k', 'd']

        core_to_vertex_map = {}
        for key, item in kd_result.groupby(['k', 'd'], as_index=False):
            assert len(item['vertex'].unique()) == item.shape[0]
            core_to_vertex_map[(int(key[0]), int(key[1]))] = list(
                map(int, item['vertex'].unique()))
            for v in core_to_vertex_map[(int(key[0]), int(key[1]))]:
                assert type(v) == int
                assert v in neighbor
        distinct_core_numbers = list(core_to_vertex_map.keys())

    
    # implementation of exp_9a
    result = {}
    
    if(args.sir_prev_9a):
        distinct_core_numbers.sort(reverse=True)
    elif(args.sir_prev_kd):
        distinct_core_numbers.sort(key=lambda x: (x[0], x[1]), reverse=True)
        core_to_vertex_map_new = {}
        distinct_core_numbers_new = []
        for c in distinct_core_numbers:
            if c[0] not in distinct_core_numbers_new:
                distinct_core_numbers_new.append(c[0])
                core_to_vertex_map_new[c[0]] = core_to_vertex_map[c]
        distinct_core_numbers = distinct_core_numbers_new
        core_to_vertex_map = core_to_vertex_map_new
    elif(args.sir_prev_dk):
        distinct_core_numbers.sort(key=lambda x: (x[1], x[0]), reverse=True)
    else:
        raise NotImplementedError()

    # print(distinct_core_numbers)
    # quit()    

    pool = multiprocessing.Pool(processes=num_rounds)
    result_all_run = pool.map(
        worker, [(distinct_core_numbers, core_to_vertex_map, neighbor)]*num_rounds)

    entry['result'] = {}
    for i, result_single_run in enumerate(result_all_run):
        entry['result'][i] = result_single_run
    entry['max propagation time'] = args.max_propagation_time
    entry['seed size'] = args.seed_size
    entry['num delete'] = args.num_delete

    result = pd.DataFrame()
    result = result.append(entry, ignore_index=True)
    result.to_csv('../output/propagation_result_prev_9a.csv', header=False,
                    index=False, mode='a')

    if(True):
        print("\n")
        print(
            ", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))
