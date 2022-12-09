import pickle
import random
import argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="enron")
parser.add_argument("-n", "--num_delete", type=int,
                    default=-1, help="how many vertices are deleted")
args = parser.parse_args()


def get_nodes(dataset, algo, num_delete):
    index = 1
    neighbor_data_filename = "sirdata_naheed_vai/" + algo + "_" + \
        dataset + "_h" + str(index) + "_" + str(num_delete) + ".csv"

    # get nodes information
    nodes = {}
    with open(neighbor_data_filename) as file:
        for line in file:
            vs = line.strip().split(",")
            assert int(vs[0]) not in nodes
            nodes[int(vs[0])] = True
            # neighbor[int(vs[0])] = list(map(int, vs[1:]))
            # for u in neighbor[int(vs[0])]:
            #     assert type(u) == int

    return nodes


dataset = args.dataset
num_delete = args.num_delete
nodes_graph = get_nodes(dataset, "graph_core", num_delete)
nodes_nbr = get_nodes(dataset, "naive_nbr", num_delete)
nodes_degree = get_nodes(dataset, "naive_degree", num_delete)


common_nodes = set(nodes_nbr.keys()).intersection(
    set(nodes_degree.keys()).intersection(
        set(nodes_graph.keys())
    )
)


def get_core(filename):
    # get core information
    core_to_vertex_map = {}
    distinct_cores = []
    with open(filename) as file:
        for line in file:
            vs = line.strip().split(",")
            vs = list(map(int, vs))
            if(vs[1] not in core_to_vertex_map):
                distinct_cores.append(vs[1])
                core_to_vertex_map[vs[1]] = [vs[0]]
            else:
                core_to_vertex_map[vs[1]].append(vs[0])
    # print(core_to_vertex_map)

    distinct_cores.sort(reverse=True)
    return core_to_vertex_map, distinct_cores


core_nbr, distinct_core_nbr = get_core(
    "sirdata_naheed_vai/core_naive_nbr_" + dataset + "_h0_" + str(num_delete) + ".csv")


potential_seeds = []
for c in distinct_core_nbr:
    # apply shuffling
    shuffled_nodes = core_nbr[c].copy()
    random.shuffle(shuffled_nodes)
    for node in shuffled_nodes:
        if(node in common_nodes):
            potential_seeds.append(node)


assert len(potential_seeds) == len(common_nodes)


assert potential_seeds[0] not in core_nbr[distinct_core_nbr[0]]

filename = "sirdata/potential_seeds_" + dataset + ".pkl"
with open(filename, "wb") as handle:
    pickle.dump(potential_seeds, handle)
