import shutil
from ast import literal_eval
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.ticker as ticker
import argparse
import itertools
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="dblp")
args = parser.parse_args()

output_folder = '../output/'

fontsize = 28
labelsize = 20
save = True
    

dataset = args.dataset
processed_filename = "../output/processed_propagation_result_9b_" + dataset + ".csv"
time_step_list = [100]
if(not os.path.isfile(processed_filename)):

    src = '../output/'
    cols = ['algo', 'dataset', 'exp', 'intervention_results', 'max propagation time',
            'num delete', 'p', 'result', 'seed size', 'timestep_results']
    # cols = ['dataset', 'p', 'algo', 'exp', 'result', 'timestep_results', 'intervention_results', 'num delete']
    # cols = ['dataset', 'p', 'algo', 'exp', 'result', 'timestep_results', 'intervention_results']

    df_nbr = pd.read_csv(
        src + "propagation_result_recursive_delinner_" + dataset + "_naive_nbr3.csv", header=None)
    df_nbr.columns = cols

    df_degree = pd.read_csv(
        src + "propagation_result_recursive_delinner_" + dataset + "_naive_degree3.csv", header=None)
    df_degree.columns = cols

    df_graph_core = pd.read_csv(
        src + "propagation_result_recursive_delinner_" + dataset + "_graph_core3.csv", header=None)
    df_graph_core.columns = cols

    df = pd.concat([df_degree, df_nbr, df_graph_core])
    df.head(n=10)


    

    # plt.style.use('grayscale')
    # from matplotlib.ticker import MaxNLocator
    # sns.set(rc={'figure.figsize': (7, 4)})
    # plt.rcParams['figure.figsize'] = (7,5)

    lw = 3
    output_folder = '../fig/'
    topk = 5

    ignore_datasets = ['bin_1', 'bin_2', 'bin_4',
                    'bin_5', 'congress', 'contact']
    group_list = ['dataset', 'p', 'algo', 'num delete']
    goodname_algo = {
        'graph_core': 'clique',
        'naive_nbr': 'nbr',
        'naive_degree': 'degree'
    }
    order = [goodname_algo[a]
            for a in ['naive_nbr', 'graph_core', 'naive_degree']]

    df_plot = None
    for key, item in df[(df['intervention_results'].notnull())].groupby(group_list, as_index=False):

        print(key)
        # print(item)
        # continue

        item['algo'] = item['algo'].replace(goodname_algo)
        assert len(item['algo'].unique()) == 1

        result = literal_eval(item['intervention_results'].iloc[0])

        result_df = pd.DataFrame()
        result_tuple = [('H'+str(hypergraph_id), k, iteration, time_step + 1, node_id, time_step_results, key[3])
                        for hypergraph_id in sorted(list(result.keys()))
                        for k in result[hypergraph_id]
                        for iteration, all_iteration_results in enumerate(result[hypergraph_id][k])
                        for node_id, single_iteration_results in enumerate(all_iteration_results)
                        for time_step, time_step_results in enumerate(single_iteration_results[2])]

        result_df = result_df.append(pd.DataFrame(result_tuple, columns=[
            'hypergraph', 'core number', 'iteration', 'time_step', 'seed id', 'infected', 'num delete']), ignore_index=False)

        # print(meandf)
        merged_df = pd.merge(result_df[result_df['hypergraph'] == "H0"], result_df[result_df['hypergraph'] == "H1"],
                             how="inner", on=["core number", "iteration", "time_step", 'seed id', 'num delete'])
        merged_df['infected difference'] = merged_df.apply(
            lambda x: x['infected_x'] - x['infected_y'], axis=1)
        merged_df['Decomposition'] = goodname_algo[key[2]]

        merged_df.drop(['infected_x', 'infected_y', 'hypergraph_x',
                        'hypergraph_y'], axis=1, inplace=True)

        # print(merged_df)

        if(df_plot is None):
            df_plot = merged_df.copy()
        else:
            df_plot = df_plot.append(merged_df, ignore_index=True)

    df_plot.to_csv(processed_filename, header=True, index=False, mode='a')
else:
    df_plot = pd.read_csv(processed_filename)

for time_step in time_step_list:

    final_legend_dic = {
        'nbr': 'Nbr',
        'degree': 'Degree',
        'clique': 'Clique',
        'kd': '(k, d)',
        "dk": '(d, k)'
    }

    hatch_dict = {
        'Nbr': "O",
        '(k, d)': '.',
        'Clique': 'x',
        'Degree': '*'
    }
    include_algos = ['Clique', 'Degree', 'Nbr']

    # df_plot['num delete'] = df_plot.apply(lambda x: 1000 if x['num delete'] == -1 else x['num_delete'], axis=1)
    df_plot['Decomposition'] = df_plot['Decomposition'].replace(final_legend_dic)

    #print(df_plot['num delete'].dtype)
    df_plot['num delete'] = df_plot['num delete'].astype(int)
    #print(df_plot['num delete'].dtype)

    # assert len(df_plot['core number sorted'].unique()) == 1
    sns.set_style("dark", {'axes.grid' : True})
    plt.style.use('grayscale')
    fig, ax = plt.subplots(figsize=(8, 4))
    bar = sns.barplot(y='infected difference', x='num delete',
                hue='Decomposition', hue_order=include_algos, data=df_plot[df_plot['time_step'] == time_step], color='k')
    plt.xlabel('#deleted nodes from inner core', fontsize=fontsize - 2)
    plt.ylabel("Decrease in avg.\nspread per seed", fontsize=fontsize - 2)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x+1) + 'K'))
    #ax.xaxis.set_major_formatter(ticker.EngFormatter())
    plt.xticks(fontsize=fontsize - 2)
    plt.yticks([0, 200, 400, 600, 800], fontsize=fontsize - 2)
    # plt.grid(axis="y")

    # h = itertools.cycle([hatch_dict[i] for i in include_algos])
    # for i,thisbar in enumerate(bar.patches):
    #     if i%len(include_algos)==0:
    #         hatch = next(h)
    #     thisbar.set_hatch(hatch)

    # plt.legend(loc='best', fontsize=fontsize-4)
    plt.legend(loc='upper center', bbox_to_anchor=(
        0.5, 1.3), ncol=4, fancybox=False, shadow=True, fontsize=labelsize-2, columnspacing=0.8)
    # plt.title(dataset + ", #delete: " + str(num_delete) +", time:" + str(time_step), fontsize=fontsize - 4)
    # plt.title(dataset + ", " + "time: " + str(time_step), fontsize=fontsize - 4)
    # leg = plt.legend(loc="best", fontsize=fontsize-6,
    #                 frameon=False, bbox_to_anchor=(0.44, 0.52))
    # for legobj in leg.legendHandles:
    #     legobj.set_linewidth(4.0)

    plt.tight_layout()
    filename = "../fig/" + dataset + "_9c_" + str(time_step) + ".pdf"
    print(filename)
    # print(output_folder)
    if(save):
        plt.savefig(filename)
        plt.show()
        pass
    else:
        # print(filename)
        # plt.show()
        pass
        # break
    plt.clf()
