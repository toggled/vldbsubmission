import sys
import os
sys.path.append("HyperNetX")
import hypernetx as hnx
import utils
import pandas as pd


verbose = True

scenes = {
    0: ('FN', 'TH'),
    1: ('TH', 'JV'),
    2: ('BM', 'FN', 'JA'),
    3: ('JV', 'JU', 'CH', 'BM'),
    4: ('JU', 'CH', 'BR', 'CN', 'CC', 'JV', 'BM'),
    5: ('TH', 'GP'),
    6: ('GP', 'MP'),
    7: ('MA', 'GP')
}

H = hnx.Hypergraph(scenes)



nodes = list(H.nodes)
num_nodes = len(nodes)
core = {}
bucket = {}  # mapping of number of neighbors, say n, to nodes with exactly n neighbors

# auxiliary data, that improves efficiency
# node_to_neighbors = {} # Not necessary, I guess
node_to_degree= {}  # inverse bucket
# removed_nodes = []


# Initial bucket fill-up
max_degree = -1
for node in nodes:

    degree = H.degree(node)
    if(degree > max_degree):
        max_degree = degree
    node_to_degree[node] = degree
    # print(node, neighbors)
    if degree not in bucket:
        bucket[degree] = [node]
    else:
        bucket[degree].append(node)

entry = {}
# add more keys here as per necessity
entry['bucket size'] = len(bucket)
entry['bucket update'] = 0

print("\n---------- Initial neighbors -------")
for node in H.nodes:
    print(node, H.neighbors(node))
print()

print("\n---------- Initial bucket -------")
print(bucket)
print()

for k in range(1, max_degree + 1):

    # TODO Discuss with Naheed vai about this case
    # if(k not in bucket):
    #     continue

    # Inner while loop
    # assert k in bucket
    while len(bucket.get(k, [])) != 0:
        v = bucket[k].pop(0)  # get first element in the
        print("k:", k, "node:", v)
        core[v] = k
        temp_nodes = nodes[:]  # (Naheed:) Make a copy of V. list(H.nodes) is not the same as V, because list(H\{
        # v}.nodes) is not the same as V\{v}
        temp_nodes.remove(v)  # V' <- V \ {v}

        # H_temp = H.restrict_to_nodes(temp_nodes) # (Naheed:) This is not what we want.
        H_temp = utils.strong_subgraph(H, temp_nodes)

        # enumerating over all neighbors of v
        # for u in H.neighbors(v):
        # for u in node_to_neighbors[v]:
        for u in utils.get_nbrs(H, v):  # (Naheed:) H.neighbors(v) does not return [] when v \notin H. Wrote a wrapper.

            print(node_to_degree)

            print("Considering neighbor", u)
            # neighbors_u = list(H_temp.neighbors(u))  # (Naheed: ) Unnecessary
            # len_neighbors_u = len(neighbors_u)  # (Naheed:) Buggy. Hypernetx neighbors() function may return None
            # print(u, 'has neighbors on sub-hypergraph:', neighbors_u)

            degree_u = utils.get_degree(H_temp, u) # Naheed vai: please check correctness of this implementation

            max_value = max(degree_u, k)
            print("max core between", k, 'and', degree_u, "is ", max_value)
            print("The location of", u, "is updated from", node_to_degree[u], "to", max_value)
            bucket[node_to_degree[u]].remove(u)

            # Move u to new location in bucket
            if max_value not in bucket:
                # TODO does code reach here?
                bucket[max_value] = [u]
            else:
                bucket[max_value].append(u)
                entry['bucket update'] += 1

            # update new location of u
            node_to_degree[u] = max_value

            print("-------- Updated bucket ---------")
            print(bucket)
            print()

        nodes = temp_nodes
        H = H_temp

print("\n\nOutput")
print(core)



# result dump
result = pd.DataFrame()
result = result.append(entry, ignore_index=True)
os.system("mkdir -p output")
result.to_csv('output/result.csv', header=False,
                        index=False, mode='a')


if(verbose):
    print("\n\nColumns in pandas") 
    print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))


