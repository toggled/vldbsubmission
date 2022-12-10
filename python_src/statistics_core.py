import numpy as np
import pandas as pd
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


# neighbor = get_neighbor(dataset)

top_core_info = {}
top_core_order = {}
top = 10
for algo in ['graph_core', 'naive_nbr', 'naive_degree']:
    core, distinct_core = get_core(
        "sirdata_naheed_vai/core_" + algo + "_" + dataset + "_h0_-1.csv")

    # get top k core info
    top_core_info[algo] = {}
    for i in range(min(top, len(distinct_core))):
        top_core_info[algo][distinct_core[i]] = core[distinct_core[i]]
    top_core_order[algo] = distinct_core[:min(top, len(distinct_core))]

print()
print("="*50)
print(dataset)
print("="*50)
for algo1, algo2 in [('naive_nbr', 'graph_core'), ('naive_nbr', 'naive_degree'), ('naive_degree', 'graph_core')]:
    result = []
    for i in range(min(len(top_core_info[algo1]), len(top_core_info[algo2]))):

        result.append((len(top_core_info[algo1][top_core_order[algo1][i]]),
                       top_core_order[algo1][i],
                       len(top_core_info[algo2]
                           [top_core_order[algo2][i]]),
                       top_core_order[algo2][i],
                       len(list(set(top_core_info[algo1][top_core_order[algo1][i]]) & set(top_core_info[algo2][top_core_order[algo2][i]])))))

    # print(result)
    result = pd.DataFrame(result, columns=['#nodes (' + algo1 + ")", 'core number (' + algo1 + ")",
                                           '#nodes (' + algo2 + ")", 'core number (' + algo2 + ")", "#common nodes"])
    result['cumulative sum'] = result['#common nodes'].cumsum()
    print()
    print("Comparison between: ", algo1, "and", algo2)
    print()
    print(result)
    print()
    result.to_csv('statistics/' + dataset + '_' + algo1 + "_" + algo2 + '.csv', header=True,
                  index=False)

print("\n"*4)

# from pprint import pprint
# print(dataset)
# pprint(result)
