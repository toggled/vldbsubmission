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

# # Visualise hypergraph for verification purposes
# hnx.drawing.draw(H,with_edge_labels = False, layout_kwargs = {'seed': 39})
# plt.savefig('test.pdf')



nodes = list(H.nodes)
num_nodes = len(nodes)
core = {}
bucket = {}  # mapping of number of neighbors, say n, to nodes with exactly n neighbors

# auxiliary data, that improves efficiency
# node_to_neighbors = {} # Not necessary, I guess
node_to_num_neighbors = {}  # inverse bucket
# removed_nodes = []


# Initial bucket fill-up
for node in nodes:
    neighbors = list(H.neighbors(node))
    len_neighbors = len(neighbors)  # this computation can be repeated
    # node_to_neighbors[node] = neighbors
    node_to_num_neighbors[node] = len_neighbors
    # print(node, neighbors)
    if len_neighbors not in bucket:
        bucket[len_neighbors] = [node]
    else:
        bucket[len_neighbors].append(node)


print("\n---------- Initial neighbors -------")
for node in H.nodes:
    print(node, H.neighbors(node))
print()

print("\n---------- Initial bucket -------")
print(bucket)
print()

for k in range(1, num_nodes + 1):

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
        H_temp = utils.strong_subgraph(H, temp_nodes) # Store.... + executation time.. 

        # enumerating over all neighbors of v
        # for u in H.neighbors(v):
        # for u in node_to_neighbors[v]:
        for u in utils.get_nbrs(H, v):  # (Naheed:) H.neighbors(v) does not return [] when v \notin H. Wrote a wrapper.

            print(node_to_num_neighbors)

            print("Considering neighbor", u)
            # neighbors_u = list(H_temp.neighbors(u))  # (Naheed: ) Unnecessary
            # len_neighbors_u = len(neighbors_u)  # (Naheed:) Buggy. Hypernetx neighbors() function may return None
            # print(u, 'has neighbors on sub-hypergraph:', neighbors_u)

            len_neighbors_u = utils.get_number_of_nbrs(H_temp, u)  # To avoid bug, we use our wrapper function.
            # How many times is neighborhood computation done? and executation time...

            max_value = max(len_neighbors_u, k)
            print("max core between", k, 'and', len_neighbors_u, "is ", max_value)
            print("The location of", u, "is updated from", node_to_num_neighbors[u], "to", max_value)
            bucket[node_to_num_neighbors[u]].remove(u)

            # Move u to new location in bucket
            if max_value not in bucket:
                # TODO does code reach here?
                bucket[max_value] = [u]
            else:
                bucket[max_value].append(u)
                entry['bucket update'] += 1

            # How many times is bucket updated + executation time??? Store...

            # update new location of u
            node_to_num_neighbors[u] = max_value

            print("-------- Updated bucket ---------")
            print(bucket)
            print()

        nodes = temp_nodes
        H = H_temp

print("\n\nOutput")
print(core) # Store...



# result dump
result = pd.DataFrame()
result = result.append(entry, ignore_index=True)
os.system("mkdir -p output")
result.to_csv('output/result.csv', header=False,
                        index=False, mode='a')


if(verbose):
    print("\n\nColumns in pandas") 
    print(", ".join(["\'" + column + "\'" for column in result.columns.tolist()]))


