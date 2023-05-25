import sys
import matplotlib

matplotlib.use('TkAgg')
sys.path.append("HyperNetX")
import hypernetx as hnx
import math
from copy import deepcopy
from hgDecompose import utils

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

glb = math.inf
llb = {}
gub = -math.inf
lub = {}
# flag = {}
# Initial bucket fill-up
for node in nodes:
    nbrs = H.neighbors(node)
    # neighbors = list(H.neighbors(node))
    len_neighbors = len(nbrs)  # this computation can be repeated
    glb = min(glb, len_neighbors)
    gub = max(gub, len_neighbors)
    # LB2 computation
    for u in nbrs:
        llb[node] = min(llb.get(node, len_neighbors), len(H.neighbors(u)) - 1)
    # node_to_neighbors[node] = neighbors
    node_to_num_neighbors[node] = len_neighbors
    # print(node, neighbors)
    if len_neighbors not in bucket:
        bucket[len_neighbors] = [node]
    else:
        bucket[len_neighbors].append(node)
    # flag[node] = False

print("\n---------- Initial neighbors -------")
for node in H.nodes:
    print(node, H.neighbors(node))
print()

print("\n---------- Initial bucket -------")
print(bucket)
print()

# Compute Local upper bounds
copy_bucket = deepcopy(bucket)
inv_bucket = deepcopy(node_to_num_neighbors)

for k in range(glb, gub):
    while len(copy_bucket.get(k, [])) != 0:
        v = copy_bucket[k].pop(0)
        lub[v] = k
        for u in utils.get_nbrs(H, v):
            if u not in lub:
                max_value = max(inv_bucket[u] - 1, k)
                copy_bucket[inv_bucket[u]].remove(u)
                copy_bucket[max_value].append(u)
                inv_bucket[u] = max_value

print('local upper bound: ')
print(sorted(lub.items()))
print('local lower bound: ')
print(sorted(llb.items()))


# Interval generator function (s is a parameter)
def generate_intervals(llb, lub, s=1):
    min_llb = min([llb[u] for u in llb])
    ub_set = set([lub[u] for u in lub]).union([min_llb - 1])
    sorted_ub_set = sorted(ub_set, reverse=True)
    i = s
    while i < len(ub_set):
        yield sorted_ub_set[i] + 1, sorted_ub_set[i - s]
        if i+s < len(ub_set):
            i += s
        else:
            if i != len(ub_set) - 1:
                yield sorted_ub_set[-1] + 1, sorted_ub_set[i]
            i += s


gen = generate_intervals(llb, lub, s=3)
final_bucket = {}
setlb = {}
inv_bucket = {}

for lower, upper in gen:
    print("Inverval [%d,%d]"%(lower, upper))
    V_kmin = [u for u in nodes if lub[u] >= lower]
    for u in V_kmin:
        if u in core:
            max_val = max(lower-1, llb[u], core[u])
        else:
            max_val = max(lower-1, llb[u])
        final_bucket[max_val] = final_bucket.get(max_val, [])+[u]
        inv_bucket[u] = max_val
        setlb[u] = True
    H_kmin = utils.strong_subgraph(H, V_kmin)

    for k in range(lower-1, upper+1):
        while len(final_bucket.get(k, [])) != 0:
            v = final_bucket[k].pop(0)
            # core[v] = k
            if setlb[v]:
                num_nbrs_v = utils.get_number_of_nbrs(H_kmin,v)
                final_bucket[num_nbrs_v] = final_bucket.get(num_nbrs_v,[])+[v]
                inv_bucket[v] = num_nbrs_v
                setlb[v] = False
            else:
                if k >= lower:
                    core[v] = k
                    setlb[v] = True
                _temp_nodes = V_kmin[:]  # Make a copy of nodes
                _temp_nodes.remove(v)  # V' <- V \ {v}
                _temp_Hkmin = utils.strong_subgraph(H_kmin, _temp_nodes)
                for u in utils.get_nbrs(H_kmin, v):
                    len_neighbors_u = utils.get_number_of_nbrs(_temp_Hkmin, u)
                    max_value = max(len_neighbors_u, k)
                    if max_value != inv_bucket[u]:
                        final_bucket[max_value] = final_bucket.get(max_value,[])+[u]
                        final_bucket[inv_bucket[u]].remove(u)
                        inv_bucket[u] = max_value
                V_kmin = _temp_nodes
                H_kmin = _temp_Hkmin

print("\n\nOutput")
print(core)
