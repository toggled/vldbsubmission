from hgDecompose.utils import get_hg, writeHypergraphHg, writeHypergraph, get_random_hg, save_dictascsv
from hgDecompose.Hypergraph import Hypergraph
import random
import os
import time
import argparse
import pickle
import copy
from hgDecompose.optimizedhgDecompose import HGDecompose
from matplotlib import pyplot as plt
import networkx as nx
# random.seed(10)
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="default")
parser.add_argument("-a", "--algo", type=str, default="naive_nbr")
parser.add_argument("-n", "--num_delete", type=int,
                    default=-1, help="how many vertices are deleted")
parser.add_argument(
    "-l", "--level", help="how many times innermost core is deleted", default=1, type=int)
args = parser.parse_args()


def construct_pickle():

    def get_core_neighbor_data(index=0):
        core_data_filename = "sirdata_backup/core_" + \
            args.algo + "_" + args.dataset + "_h" + str(index) + ".csv"

        # get core information
        core_base = {}
        with open(core_data_filename) as file:
            for line in file:
                vs = line.strip().split(",")
                assert int(vs[0]) not in core_base
                core_base[int(vs[0])] = int(vs[1])

        neighbor_data_filename = "sirdata_backup/log_" + \
            args.dataset + "_h" + str(index) + "_" + \
            str(args.num_delete) + ".csv"

        # get neighborhood information
        neighbor = {}
        with open(neighbor_data_filename) as file:
            for line in file:
                vs = line.strip().split(",")
                assert int(vs[0]) not in neighbor
                neighbor[int(vs[0])] = list(map(int, vs[1:]))
                for u in neighbor[int(vs[0])]:
                    assert type(u) == int
                assert int(vs[0]) in core_base, str(int(vs[0])) + " " + str(index)

        return core_base, neighbor

    output = {}

    for index in [0, 1]:
        output[index] = {}
        output[index]['H'] = None
        output[index]['core'], output[index]['neighbor'] = get_core_neighbor_data(
            index=index)

    # store as pickle
    pathstring = "sirdata/"
    os.system('mkdir -p ' + pathstring)
    os.system("mkdir -p tests/tmp")
    with open(os.path.join(pathstring, args.dataset + '_' + args.algo + "_" + str(args.num_delete) + '.pkl'), 'wb') as handle:
        pickle.dump(output, handle, protocol=4)


construct_pickle()
