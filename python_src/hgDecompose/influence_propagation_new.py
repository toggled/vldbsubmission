import random
from tqdm import tqdm


def propagate_for_all_vertices(neighbor, core, num_vertex_per_core=100, top_k=5,  p=0.5, num_iterations=100, original_n=None, verbose=True):

    result = {}  # Entry is a core number. value is a list of percentages of the infected population for all vertices with the same core number

    core_to_vertex_map = {}
    distinct_core_numbers = []
    for v in core:
        if(core[v] not in core_to_vertex_map):
            core_to_vertex_map[core[v]] = [v]
            distinct_core_numbers.append(core[v])
        else:
            core_to_vertex_map[core[v]].append(v)

    distinct_core_numbers.sort(reverse=True)

    for core_number in tqdm(distinct_core_numbers[:top_k]):
        # check size
        v_sampled = None
        if(len(core_to_vertex_map[core_number]) < num_vertex_per_core):
            v_sampled = core_to_vertex_map[core_number]
        else:
            v_sampled = random.choices(
                core_to_vertex_map[core_number], k=num_vertex_per_core)

        for v in tqdm(v_sampled):
            if(core_number not in result):
                result[core_number] = [bfs_bounded(
                    neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose)]
            else:
                result[core_number].append(bfs_bounded(
                    neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose))

    return result


def propagate_for_all_vertices_for_kd(neighbor, kd_core, num_vertex_per_core=100, top_k=5,  p=0.5, num_iterations=10, original_n=None, verbose=True):

    result = {}  # Entry is a core number. value is a list of percentages of the infected population for all vertices with the same core number

    distinct_core_numbers = sorted(
        kd_core.keys(), key=lambda tup: (tup[0], tup[1]), reverse=True)

    # print(distinct_core_numbers)
    visited_core_number = []
    for core_number in tqdm(distinct_core_numbers):
        # allow top_k core
        if(core_number[0] in visited_core_number):
            continue
        if(len(visited_core_number) >= top_k):
            break

        visited_core_number.append(core_number[0])
        # print(core_number)

        # check size
        v_sampled = None
        if(len(kd_core[core_number]) < num_vertex_per_core):
            v_sampled = kd_core[core_number]
        else:
            v_sampled = random.choices(
                kd_core[core_number], k=num_vertex_per_core)
        # print(v_sampled)
        for v in tqdm(v_sampled):
            if(core_number not in result):
                result[core_number] = [bfs_bounded(
                    neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose)]
            else:
                result[core_number].append(bfs_bounded(
                    neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose))
            # print(component_sz(v))

    return result


def bfs_bounded(neighbor, starting_vertex, p=0.3, num_iterations=10, original_n=None, verbose=True):
    visited = {}
    for v in neighbor.keys():
        visited[v] = False
    results_each_time_step = []

    queue = []
    queue.append(starting_vertex)
    visited[starting_vertex] = True
    num_infected = 1
    i = 0
    while queue and i < num_iterations:
        i += 1
        v = queue.pop(0)
        for u in neighbor[v]:
            if(not visited[u]):
                if(random.random() <= p):
                    queue.append(u)
                    visited[v] = True
                    num_infected += 1

        results_each_time_step.append(num_infected)

    return num_infected, len(neighbor[starting_vertex]), results_each_time_step


def run_intervention_exp2(name, original_n, p=0.5, top_k=5, verbose=False):
    import os
    import pickle
    import numpy as np
    pkl_path = 'sirdata/'
    path = pkl_path+name+'.pkl'
    with open(os.path.join(path), 'rb') as handle:
        data = pickle.load(handle)
        print("loaded ", path)
    result = {}
    for k in data:
        print('Core deletion#: ', k)
        result[k] = {}
        temp_core = data[k]['core']
        neighbor = data[k]['neighbor']
        # check_connectivity(H)
        # result[k] = propagate_for_all_vertices(H, temp_core, p = p, original_n = original_n, verbose=verbose)
        core_to_vertex_map = {}
        distinct_core_numbers = []
        for v in temp_core:
            if(temp_core[v] not in core_to_vertex_map):
                core_to_vertex_map[temp_core[v]] = [v]
                distinct_core_numbers.append(temp_core[v])
            else:
                core_to_vertex_map[temp_core[v]].append(v)

        distinct_core_numbers.sort(reverse=True)
        for core_number in distinct_core_numbers[:top_k]:
            print('core: ', core_number)
            # result[k][core_number] = {}
            if(len(core_to_vertex_map[core_number]) < 100):
                v_sampled = core_to_vertex_map[core_number]
            else:
                v_sampled = random.choices(
                    core_to_vertex_map[core_number], k=100)

            num_iterations = 100
            result_all_run = []
            for _ in range(5):
                result_single_run = []
                for v in v_sampled:
                    result_single_run.append(bfs_bounded(
                        neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose)[0])
                result_all_run.append(result_single_run)
            result[k][core_number] = list(
                np.array(result_all_run).mean(axis=0))
    return result
