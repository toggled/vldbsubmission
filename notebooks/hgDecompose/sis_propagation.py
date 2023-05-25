import random
import numpy as np


def propagateSIS_for_all_vertices(H, core, p = 0.5, gamma = 0.01, num_iterations = 2, verbose=True):

    result = {} # Entry is a core number 
    # value is a list of percentages of the infected population for all vertices with the same core number
    for key in core:
        if (verbose):
            print("\n====== core-number: =====>", core[key], ' vertex: ',key)
        if(core[key] not in result):
            result[core[key]] = [propagateSIS(H, starting_vertex=key, p = p, num_iterations = num_iterations, verbose = verbose)]
        else:
            result[core[key]].append(propagateSIS(H, starting_vertex=key, p = p, num_iterations = num_iterations, verbose = verbose))

    return result

def propagateSIS(H, starting_vertex, p = 0.5, gamma = 0.01, num_iterations = 2, verbose=True):
    """
    SIS epidemic spread
    """
    random.seed(10)
    suscepted = H.nodes()
    suscepted.remove(starting_vertex)
    infected = [starting_vertex]
    recovered = []

    for i in range(num_iterations):
        if(verbose):
            print('\n\n\nIteration:', i)
            print("infected:", infected)
            print("recovered:", recovered)
            print("suscepted:", suscepted)
            print()
        
        if(len(infected) == 0):
            if(verbose):
                print("No more propagation is possible")
            break
        
        new_infected = []
        new_recovered = []    
        for v in infected:
            if(verbose):
                print("\nPorpagating for", v)
            for u in H.neighbors(v):
                if(u in suscepted): # a neighbour of infected v: susceptible -> infected with prob p.
                    if(random.random() <= p):
                        if(verbose):
                            print(v, "->", u)
                        new_infected.append(u)
                        suscepted.remove(u)
                    else:
                        if(verbose):
                            print(v, "->", u, "not propagated")
                # else:
                #     if(verbose):
                #         print(u, "is already either infected or recovered")
            new_recovered.append(v)

        infected += new_infected
        recovered += new_recovered
        for v in new_recovered:
            infected.remove(v)

        """ Re-infect a recovered person with gamma probability """
        new_infected = []
        for v in recovered[:]:
            if(random.random() <= gamma):
                if (verbose):
                    print(v, 'is reinfected')
                new_infected.append(v)
                recovered.remove(v)
        infected += new_infected

    return 1 - float(len(suscepted) / H.get_N())
