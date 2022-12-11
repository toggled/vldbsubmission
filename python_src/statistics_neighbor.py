import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="enron")
# parser.add_argument("-a", "--algo", type=str, default="naive_nbr")
# parser.add_argument("-n", "--num_delete", type=int,
#                     default=-1, help="how many vertices are deleted")
# parser.add_argument(
#     "-l", "--level", help="how many times innermost core is deleted", default=1, type=int)
args = parser.parse_args()


dataset = args.dataset

result = {}
result_sum ={}

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


def get_neighbor(dataset):
    neighbor_data_filename = "sirdata_naheed_vai/naive_nbr_" + dataset + "_h0_-1.csv"
    neighbor = {}
    with open(neighbor_data_filename) as file:
        for line in file:
            vs = line.strip().split(",")
            assert int(vs[0]) not in neighbor
            neighbor[int(vs[0])] = list(map(int, vs[1:]))
    return neighbor


neighbor = get_neighbor(dataset)
from scipy import stats as st
for algo in ['graph_core', 'naive_nbr', 'naive_degree']:
    core, distinct_core = get_core(
        "sirdata_naheed_vai/core_" + algo + "_" + dataset + "_h0_-1.csv")

    print("innermost core number", algo, distinct_core[0], "innermost core length", len(core[distinct_core[0]]))
    neighbor_len = []
    for node in core[distinct_core[0]]:
        neighbor_len.append(len(neighbor[node]))
    assert len(neighbor_len) == len(core[distinct_core[0]])


    avg_neighbor_len = np.array(neighbor_len).mean()
    # avg_neighbor_len = st.mode(neighbor_len)[0]
    result[algo] = avg_neighbor_len
    result_sum[algo] = avg_neighbor_len * len(core[distinct_core[0]])

from pprint import pprint
print(dataset)
pprint(result)

print()

pprint(result_sum)
