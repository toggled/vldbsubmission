import multiprocessing
import argparse
import pandas as pd
from tqdm import tqdm
from hgDecompose.influence_propagation_new import bfs_bounded
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="enron")
parser.add_argument("-a", "--algo", type=str, default="naive_nbr")
parser.add_argument("-n", "--num_delete", type=int,
                    default=-1, help="how many vertices are deleted")
parser.add_argument(
    "-p", "--prob", help="parameter for Probability", default=0.3, type=float)
args = parser.parse_args()

seed_size = 100
num_rounds = 2
max_propagation_time = 100

def worker(arg_tuple):
    # print(arg_tuple)
    potential_seeds, neighbor = arg_tuple[0], arg_tuple[1]
    result_single_run = []
    for v in tqdm(potential_seeds[:seed_size]):
        assert v in neighbor, v
        result_single_run.append(bfs_bounded(
            neighbor, starting_vertex=v, p=args.prob, num_iterations=max_propagation_time, verbose=False))
    # print('I am number %d in process %d' % (procnum, getpid()))
    return result_single_run


if __name__ == '__main__':
    import os
    import pickle
    import numpy as np
    pkl_path = 'sirdata/'
    result = {}

    filename = pkl_path + "potential_seeds_" + args.dataset + ".pkl"
    potential_seeds = None
    with open(filename, "rb") as handle:
        potential_seeds = pickle.load(handle)

    for k in tqdm(range(2), disable=True):
        print('Core deletion#: ', k)
        result[k] = {}
        neighbor = {}
        neighbor_data_filename = "sirdata_naheed_vai/" + args.algo + \
            "_" + args.dataset + "_h" + \
            str(k) + "_" + str(args.num_delete) + ".csv"
        neighbor = {}
        with open(neighbor_data_filename) as file:
            for line in file:
                vs = line.strip().split(",")
                assert int(vs[0]) not in neighbor
                neighbor[int(vs[0])] = list(map(int, vs[1:]))

        for v in potential_seeds:
            assert v in neighbor, v

        # print(len(potential_seeds), len(neighbor))

        pool = multiprocessing.Pool(processes = num_rounds)
        result_all_run = pool.map(worker, [(potential_seeds, neighbor)]*num_rounds)
        # print(val)

        # result_all_run = []
        # for _ in range(3):
            
        #     result_all_run.append(result_single_run)
        # print(result_all_run)
        result[k][0] = result_all_run
    # return result

    # pool = multiprocessing.Pool(processes = 2)
    # val = pool.map(worker, range(4))
    # print(val)

    entry = {}
    entry['dataset'] = args.dataset
    entry['p'] = float(args.prob)
    entry['algo'] = args.algo
    entry['exp'] = "sir_exp3"
    entry['result'] = None
    entry['timestep_results'] = None
    entry['intervention_results'] = None
    entry['max propagation time'] = None



    entry['intervention_results'] = result
    entry['num delete'] = args.num_delete
    result = pd.DataFrame()
    result = result.append(entry, ignore_index=True)
    result.to_csv('../output/propagation_result_recursive_delinner_'+args.dataset+"_"+args.algo+'3.csv', header=False,
                  index=False, mode='a')
