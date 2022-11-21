import random
import numpy as np
from copy import deepcopy
from multiprocessing import Pool
import os
import pickle
from hgDecompose.utils import check_connectivity, component_sz, save_dict, avg_shortest_pathlen
pkl_path = 'sirdata/'


def propagate_for_all_vertices(H, core, num_vertex_per_core=100, top_k=100,  p=0.5, num_iterations=100, original_n=None, verbose=True):

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



    for core_number in distinct_core_numbers[:top_k]:
        # check size
        v_sampled = None
        if(len(core_to_vertex_map[core_number]) < num_vertex_per_core):
            v_sampled = core_to_vertex_map[core_number]
        else:
            v_sampled = random.choices(core_to_vertex_map[core_number], k=num_vertex_per_core)

        

        for v in v_sampled:
            if(core_number not in result):
                result[core_number] = [propagate2(
                    H, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose)]
            else:
                result[core_number].append(propagate2(
                    H, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose))
            # print(component_sz(v))

    # TODO: Parallelize this loop
    # core_v_list = []
    # core_numbers = []
    # for core_number in distinct_core_numbers[:top_k]:
    #     for v in random.choices(core_to_vertex_map[core_number], k=num_vertex_per_core):
    #         core_v_list.append((H,v,p,num_iterations,verbose))
    #         core_numbers.append(core_number)

    # with Pool(processes=8) as pool:
    #     pool_results = pool.map(propagate, core_v_list)
    #     # for x in pool.map(propagate, core_v_list):
    #     #     pool_results.append(x[0])
    #     pool.join()

    # for i, core_number in enumerate(core_numbers):
    #     if(core_number not in result):
    #         result[core_number] = pool_results[i][0]
    #     else:
    #         result[core_number].append(pool_results[i][0])
    #     # if(core_number not in result):
    #     #     result[core_number] = [propagate(H, starting_vertex=v, p = p, num_iterations = num_iterations, verbose = verbose)[0]]
    #     # else:
    #     #     result[core_number].append(propagate(H, starting_vertex=v, p = p, num_iterations = num_iterations, verbose = verbose)[0])
    return result

def propagate_for_all_vertices_for_kd(H, kd_core, num_vertex_per_core=100, top_k=150,  p=0.5, num_iterations=10, original_n=None, verbose=True):

    result = {}  # Entry is a core number. value is a list of percentages of the infected population for all vertices with the same core number

    distinct_core_numbers = sorted(kd_core.keys(), key=lambda tup: (tup[0], tup[1]), reverse=True)
    
    for core_number in distinct_core_numbers[:top_k]:

        # check size
        v_sampled = None
        if(len(kd_core[core_number]) < num_vertex_per_core):
            v_sampled = kd_core[core_number]
        else:
            v_sampled = random.choices(kd_core[core_number], k=num_vertex_per_core)

        for v in v_sampled:
            if(core_number not in result):
                result[core_number] = [propagate2(
                    H, starting_vertex=str(v), p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose)]
            else:
                result[core_number].append(propagate2(
                    H, starting_vertex=str(v), p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose))
            # print(component_sz(v))
   
    return result


def run_intervention_exp2(name, original_n, p=0.5, verbose=False):
    path = pkl_path+name+'.pkl'
    with open(os.path.join(path), 'rb') as handle:
        data = pickle.load(handle)
        print("loaded ", path)
    result = {}
    for k in data:
        print('Core deletion#: ', k)
        result[k] = {}
        temp_core = data[k]['core']
        H = data[k]['H']
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
        for core_number in distinct_core_numbers[:5]:
            print('core: ', core_number)
            # result[k][core_number] = {}
            if(len(core_to_vertex_map[core_number]) < 100):
                v_sampled = core_to_vertex_map[core_number]
            else:
                v_sampled = random.choices(core_to_vertex_map[core_number], k=100)

            num_iterations = 10
            result_all_run = []
            for _ in range(5):
                result_single_run = []
                for v in v_sampled:
                    result_single_run.append(propagate2(
                        H, starting_vertex=v, p=p, num_iterations=num_iterations, original_n=original_n, verbose=verbose)[0])
                result_all_run.append(result_single_run)
            result[k][core_number] = list(
                np.array(result_all_run).mean(axis=0))
    return result


def run_intervention_exp2_explain(name, verbose=False):
    """ COnnected components/reachable nodes """
    # 'temp.pkl' =>
    path = pkl_path+name+'.pkl'
    with open(os.path.join(path), 'rb') as handle:
        data = pickle.load(handle)
        print("loaded ", path)
    result = {}
    for k in data:
        # if (k!=2):
        #     continue
        print('Core deletion#: ', k)
        result[k] = {}
        temp_core = data[k]['core']
        H = data[k]['H']
        print('N: ', len(H.inc_dict))
        # check_connectivity(H)
        # continue
        core_to_vertex_map = {}
        distinct_core_numbers = []
        for v in temp_core:
            if(temp_core[v] not in core_to_vertex_map):
                core_to_vertex_map[temp_core[v]] = [v]
                distinct_core_numbers.append(temp_core[v])
            else:
                core_to_vertex_map[temp_core[v]].append(v)

        distinct_core_numbers.sort(reverse=True)

        for core_number in distinct_core_numbers[:100]:
            print('core: ', core_number)
            # result[k][core_number] = {}
            tmp = {}
            for v in random.choices(core_to_vertex_map[core_number], k=100):
                # if(core_number not in result):
                #     result[core_number] = [propagate(H, starting_vertex=v, p = p, num_iterations = 100, original_n = original_n, verbose = verbose)[0]]
                # else:
                #     result[core_number].append(propagate(H, starting_vertex=v, p = p, num_iterations = 100, original_n = original_n, verbose = verbose)[0])
                tmp[v] = component_sz(v, H)
            result[k][core_number] = np.mean(list(tmp.values()))
    save_dict(result, '../output/'+name+'_comp3.pkl')


def run_intervention_exp2_explain_splen(name, verbose=False):
    """ SHortest path length """
    path = pkl_path+name+'.pkl'
    with open(os.path.join(path), 'rb') as handle:
        data = pickle.load(handle)
        print("loaded ", path)

    result = {}

    # record samples for the innermost core
    assert 0 in data  # to check if H0 is in data
    H0_core = data[0]['core']
    constant_M = data[0]['H'].get_M()
    H0_V = list(data[0]['H'].inc_dict.keys())

    print("Max sp:", constant_M)

    for k in data:
        print("\n\n\n", k)
        # if (k!=2):
        #     continue
        result[k] = {}
        temp_core = data[k]['core']
        H = data[k]['H']
        print('N: ', len(H.inc_dict))
        # check_connectivity(H)
        # continue
        core_to_vertex_map = {}
        distinct_core_numbers = []
        for v in temp_core:
            if(temp_core[v] not in core_to_vertex_map):
                core_to_vertex_map[temp_core[v]] = [v]
                distinct_core_numbers.append(temp_core[v])
            else:
                core_to_vertex_map[temp_core[v]].append(v)

        distinct_core_numbers.sort(reverse=True)

        for core_number in distinct_core_numbers[:200]:
            # print('core: ',core_number)
            # result[k][core_number] = {}
            tmp = {}
            for v in random.choices(core_to_vertex_map[core_number], k=100):
                # conditional sampling
                # if(H0_core[v] != core_number):
                #     # print("Do not consider", v)
                #     continue

                # if(core_number not in result):
                #     result[core_number] = [propagate(H, starting_vertex=v, p = p, num_iterations = 100, original_n = original_n, verbose = verbose)[0]]
                # else:
                #     result[core_number].append(propagate(H, starting_vertex=v, p = p, num_iterations = 100, original_n = original_n, verbose = verbose)[0])
                # result[k][core_number][v] = avg_shortest_pathlen(v,H,100)
                tmp[v] = avg_shortest_pathlen(v, H, 100, H0_V, constant_M)
                if (verbose):
                    print('v ', v, ' avg SP length: ',
                          result[k][core_number][v])

            # print(tmp)
            result[k][core_number] = np.mean(list(tmp.values()))

    save_dict(result, '../output/'+name+'_sp4.pkl')


def run_intervention_exp(H, core, p=0.5, verbose=False):
    # print(core)
    # deleted_ids = [2693,2804,3865,1547,2102,2960,2537, 3446, 2120, 2673]
    max_core_number = -1
    for v in core:
        if(max_core_number < core[v]):
            max_core_number = core[v]

    # print(max_core_number)

    nodes_with_max_core = []
    for v in core:
        if(core[v] == max_core_number):
            nodes_with_max_core.append(v)

    all_nodes = H.nodes()
    result = {}
    # print(all_nodes)
    strongly_induced_eids = H.get_stronglyinduced_edgeIds(nodes_with_max_core)
    if verbose:
        print('# potential edges to delete: ', len(strongly_induced_eids))

    # Top-k edge-deletion intervention strategy.
    temp_core = {}
    for node in H.nodes():
        temp_core[node] = core[node]
        # print(eid, H.get_edge_byindex(eid), temp_H.nodes(), len(temp_H.nodes()))
    result['nill'] = propagate_for_all_vertices(
        H, temp_core, p=p, original_n=len(all_nodes), verbose=verbose)

    # Random k% edge-deletion intervention strategy TODO
    for eperc in [5, 10, 15]:
        len_dele = int(len(temp_core) * eperc/100.0)
        print('will delete ', eperc, '% = ', len_dele, ' edges')
        todelete = []
        for eid in random.choices(strongly_induced_eids, k=len_dele):
            todelete.append(eid)

        print('List : ', todelete)
        temp_H = deepcopy(H)
        for eid in todelete:
            # print('deleting ',eid)
            temp_H.del_edge(eid)
        temp_core = {}
        for node in temp_H.nodes():
            temp_core[node] = core[node]
        result['top'+str(eperc)+'%'] = propagate_for_all_vertices(temp_H,
                                                                  temp_core, p=p, original_n=len(all_nodes), verbose=verbose)

    return result


def propagate_for_random_seeds(H, core, seed_size=1000, p=0.5, num_iterations=100, verbose=False):

    # print(core)
    result = {}
    #
    for v in random.choices(H.nodes(), k=seed_size):
        # print(v)
        _, timestep_of_infection = propagate(
            H, starting_vertex=v, p=p, num_iterations=num_iterations, verbose=False)
        # print(timestep_of_infection)
        # print()
        for u in timestep_of_infection:
            if(core[u] not in result):
                result[core[u]] = [timestep_of_infection[u]]
            else:
                result[core[u]].append(timestep_of_infection[u])

    # print(result)

    return result


def propagate(H, starting_vertex, p=0.5, num_iterations=10, original_n=None, verbose=True):
    """
    Returns fraction of infected
    """
    # print('original_n: ',original_n)
    timestep_of_infection = {}
    if original_n is None:
        len_nodes = H.get_N()
    else:
        len_nodes = original_n
    for v in H.nodes():
        timestep_of_infection[v] = num_iterations + 1
    suscepted = H.nodes()
    suscepted.remove(starting_vertex)
    infected = [starting_vertex]
    timestep_of_infection[starting_vertex] = 0
    recovered = []

    for i in range(num_iterations):
        if(verbose):
            print('\n\n\nIteration:', i)
            # print("infected:", infected)
            # print("recovered:", recovered)
            # print("suscepted:", suscepted)
            print()

        if(len(infected) == 0):
            # if(verbose):
            #     print("No more propagation is possible")
            break

        new_infected = []
        new_recovered = []
        for v in infected:
            # if(verbose):
            #     print("\nPorpagating for", v)
            for u in H.neighbors(v):
                if(u in suscepted):
                    if(random.random() <= p):
                        # if(verbose):
                        #     print(v, "->", u)
                        new_infected.append(u)
                        timestep_of_infection[u] = i + 1
                        suscepted.remove(u)
                    else:
                        # if(verbose):
                        #     print(v, "->", u, "not propagated")
                        pass
                # else:
                #     if(verbose):
                #         print(u, "is already either infected or recovered")
            new_recovered.append(v)

        infected += new_infected
        recovered += new_recovered
        for v in new_recovered:
            infected.remove(v)

    return 1 - float(len(suscepted) / len_nodes), timestep_of_infection


def propagate2(H, starting_vertex, p=0.5, num_iterations=10, original_n=None, verbose=True):
    """
    Returns number of infected
    """
    # print('original_n: ',original_n)
    timestep_of_infection = {}
    if original_n is None:
        len_nodes = H.get_N()
    else:
        len_nodes = original_n
    for v in H.nodes():
        timestep_of_infection[v] = num_iterations + 1
    suscepted = H.nodes()
    suscepted.remove(starting_vertex)
    infected = [starting_vertex]
    timestep_of_infection[starting_vertex] = 0
    recovered = []

    # print(starting_vertex, len(H.neighbors(starting_vertex)))
    # quit()

    for i in range(num_iterations):
        if(verbose):
            print('\n\n\nIteration:', i)
            # print("infected:", infected)
            # print("recovered:", recovered)
            # print("suscepted:", suscepted)
            print()

        if(len(infected) == 0):
            # if(verbose):
            #     print("No more propagation is possible")
            break

        new_infected = []
        new_recovered = []
        for v in infected:
            # if(verbose):
            #     print("\nPorpagating for", v)
            for u in H.neighbors(v):
                if(u in suscepted):
                    if(random.random() <= p):
                        # if(verbose):
                        #     print(v, "->", u)
                        new_infected.append(u)
                        timestep_of_infection[u] = i + 1
                        suscepted.remove(u)
                    else:
                        # if(verbose):
                        #     print(v, "->", u, "not propagated")
                        pass
                # else:
                #     if(verbose):
                #         print(u, "is already either infected or recovered")
            new_recovered.append(v)

        infected += new_infected
        recovered += new_recovered
        for v in new_recovered:
            infected.remove(v)
    # print(len_nodes, suscepted, len(suscepted))
    # print(len_nodes - len(suscepted), timestep_of_infection)
    # return len_nodes - len(suscepted), timestep_of_infection
    return len_nodes - len(suscepted), len(H.neighbors(starting_vertex))
