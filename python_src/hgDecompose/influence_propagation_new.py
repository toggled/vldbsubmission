import random
from tqdm import tqdm


def exp_9a(neighbor, core, num_vertex_per_core=100, top_k=5,  p=0.5, num_iterations=100, original_n=None, verbose=True):

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

    # repeat
    for i in range(3):

        # sort nodes according to core number
        sorted_nodes = []
        for c in distinct_core_numbers:
            # apply shuffling
            shuffled_nodes = core_to_vertex_map[c].copy()
            random.shuffle(shuffled_nodes)
            sorted_nodes += shuffled_nodes

        # propagate
        for v in tqdm(sorted_nodes[:num_vertex_per_core]):
            if(i not in result):
                result[i] = [bfs_bounded(
                    neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose)]
            else:
                result[i].append(bfs_bounded(
                    neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose))

    return result


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
    result = {}
    if(False):
        path = pkl_path+name+'.pkl'
        with open(os.path.join(path), 'rb') as handle:
            data = pickle.load(handle)
            print("loaded ", path)
        for k in tqdm(data):
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
            for core_number in tqdm(distinct_core_numbers[:top_k]):
                # print('core: ', core_number)
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
                    for v in tqdm(v_sampled):
                        result_single_run.append(bfs_bounded(
                            neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose)[0])
                    result_all_run.append(result_single_run)
                result[k][core_number] = list(
                    np.array(result_all_run).mean(axis=0))
    else:
        filename = pkl_path + "potential_seeds_" + name.split("_")[0] + ".pkl"
        potential_seeds = None
        with open(filename, "rb") as handle:
            potential_seeds = pickle.load(handle)

        # # revise potential seeds
        # assert 1 in data
        # potential_seeds_new = []
        # for v in potential_seeds:
        #     if v in data[1]['neighbor']:
        #         potential_seeds_new.append(v)
        # # print(len(potential_seeds), len(potential_seeds_new))
        # potential_seeds = potential_seeds_new
        # # quit()

        for k in tqdm(range(2)):
            print('Core deletion#: ', k)
            result[k] = {}
            neighbor = {}
            splitted_name = name.split("_")
            neighbor_data_filename = "sirdata_naheed_vai/" + splitted_name[1] + "_" + splitted_name[2] + \
                "_" + splitted_name[0] + "_h" + \
                str(k) + "_" + splitted_name[3] + ".csv"
            neighbor = {}
            with open(neighbor_data_filename) as file:
                for line in file:
                    vs = line.strip().split(",")
                    assert int(vs[0]) not in neighbor
                    neighbor[int(vs[0])] = list(map(int, vs[1:]))

            for v in potential_seeds:
                assert v in neighbor, v

            num_iterations = 100
            result_all_run = []
            for _ in range(10):
                result_single_run = []
                for v in tqdm(potential_seeds[:100]):
                    assert v in neighbor, v
                    result_single_run.append(bfs_bounded(
                        neighbor, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose))
                result_all_run.append(result_single_run)
            result[k][0] = result_all_run
    return result
