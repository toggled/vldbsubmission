import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from ast import literal_eval
import shutil
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="dblp")
args = parser.parse_args()


dataset = args.dataset


output_folder = '../output/'

fontsize = 18
labelsize = 14




processed_filename = '../output/processed_propagation_result_prev_9a_' + dataset + '.csv'
if(not os.path.isfile(processed_filename)):
    df = pd.read_csv(output_folder + "propagation_result_prev_9a.csv", header=None)
    df.columns = ['algo', 'dataset', 'exp', 'intervention_results', 'max propagation time', 'num delete', 'p', 'result', 'seed size', 'timestep_results']
    df['seed size'] = df['seed size'].astype(int)
    df['max propagation time'] = df['max propagation time'].astype(int)
    print(df.shape)
    df.head(15)
    
    df['seed size'] = df['seed size'].astype(int)
    df['max propagation time'] = df['max propagation time'].astype(int)

    goodname_algo = {
        'naive_nbr': 'nbr',
        'naive_degree': 'degree',
        'graph_core': 'clique'
    }
    df2 = df.copy()
    group_list = ['dataset']
    for key, item in df2.groupby(group_list, as_index=False):
        processed_filename = '../output/processed_propagation_result_prev_9a_' + key + '.csv'

        if(os.path.isfile(processed_filename)):
            continue

        print("processing:", processed_filename)

        item['algo'] = item['algo'].replace(goodname_algo)
        result_df = pd.DataFrame()
        group_list2 = ['algo', 'max propagation time', 'exp', 'seed size']
        for key2, item2 in item.groupby(group_list2, as_index=False):
            print(key2)
            assert item2.shape[0] == 1
            result = literal_eval(item2['result'].iloc[0])
            result = [(k, core_number, v[0], v[1], key2[0], key2[1], key2[2])
                        for k in result 
                            for core_number in result[k] 
                                for v in result[k][core_number]]
            result_df = result_df.append(pd.DataFrame(result, columns=[
                'iteration', 'core number', 'infected', 'neighbors', 'algo', 'max propagation time', 'exp']), ignore_index=False)

        result_df = result_df.groupby(
            ['iteration', 'core number', 'algo', 'max propagation time', 'exp']).mean().reset_index()
        result_df.to_csv(processed_filename, header=True,
                            index=False, mode='a')
        print(result_df)
        


final_legend_dic = {
    'nbr': 'Nbr',
    'degree': 'Degree',
    'clique': 'Clique',
    'kd': '(k, d)',
    "dk": '(d, k)'
}

import matplotlib.ticker as ticker
save = True

processed_filename = '../output/processed_propagation_result_prev_9a_' + dataset + '.csv'
result_df = pd.read_csv(processed_filename)
result_df['algo'] = result_df['algo'].replace(final_legend_dic)
for y_var in ['infected', 'neighbors'][:1]:

    # if(dataset == "dblp"):
    #     result_df = result_df[result_df['seed size'] >= 300]
    
    sns.set_style("dark", {'axes.grid' : True})
    plt.style.use('grayscale')

    hatch_dict = {
        'Nbr': "O",
        '(k, d)': '.',
        'Clique': 'x',
        'Degree': '*'
    }

    include_algos = ['Clique', 'Degree', 'Nbr', '(k, d)']
    
    fig, ax = plt.subplots(figsize=(20, 4))
    bar = sns.barplot(x='core number', y=y_var, hue='algo', hue_order=include_algos, data=result_df[result_df['core number'] <= 13], color='k')
    plt.xlabel('Core-number of seed nodes', fontsize=fontsize-2)
    plt.ylabel("#avg infected nodes", fontsize=fontsize-2)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x/1000) + 'K'))
    # ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x/1000) + 'K'))
    plt.xticks(fontsize=fontsize-2)
    # plt.yticks([0, 2000, 4000, 6000, 8000], fontsize=fontsize-2)
    plt.legend(loc='upper center', bbox_to_anchor=(
        0.5, 1.3), ncol=4, fancybox=False, shadow=True, fontsize=labelsize-2, columnspacing=0.8)

    # h = itertools.cycle([hatch_dict[i] for i in include_algos])
    # for i,thisbar in enumerate(bar.patches):
    #     if i%len(include_algos)==0:
    #         hatch = next(h)
    #     thisbar.set_hatch(hatch)


    plt.tight_layout()
    filename = dataset + "_prev_" + y_var + ".pdf"
    print(filename)
    if(save):
        plt.savefig("../fig/" + filename)
        plt.show()
    plt.clf()
    # break