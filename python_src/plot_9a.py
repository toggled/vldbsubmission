import pandas as pd
import seaborn as sns
# %matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import os
from ast import literal_eval
import shutil
import argparse
import matplotlib.ticker as ticker
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=str, default="dblp")
args = parser.parse_args()


dataset = args.dataset

output_folder = '../output/'

fontsize = 28
labelsize = 20
# plt.style.use('grayscale')
lw = 3
save = True


processed_filename = '../output/processed_propagation_result_9a_' + dataset + '.csv'
if(not os.path.isfile(processed_filename)):

    df = pd.read_csv(output_folder + "propagation_result_9a.csv", header=None)
    df.columns = ['algo', 'dataset', 'exp', 'intervention_results',
                  'max propagation time', 'num delete', 'p', 'result', 'seed size', 'timestep_results']
    df['seed size'] = df['seed size'].astype(int)
    df['max propagation time'] = df['max propagation time'].astype(int)
    print(df.shape)
    print(df.head(15))

    goodname_algo = {
        'naive_nbr': 'nbr',
        'naive_degree': 'degree',
        'graph_core': 'clique'
    }
    df2 = df.copy()
    group_list = ['dataset']
    for key, item in df2.groupby(group_list, as_index=False):
        processed_filename = '../output/processed_propagation_result_9a_' + key + '.csv'

        if(os.path.isfile(processed_filename)):
            continue

        print("processing:", processed_filename)

        item['algo'] = item['algo'].replace(goodname_algo)
        result_df = pd.DataFrame()
        group_list2 = ['algo', 'max propagation time', 'exp', 'seed size']
        for key2, item2 in item.groupby(group_list2, as_index=False):
            assert item2.shape[0] == 1
            result = literal_eval(item2['result'].iloc[0])
            result = [(k, item2['seed size'].item(), v[0], v[1], key2[0], key2[1], key2[2])
                      for k in result for v in result[k]]
            result_df = result_df.append(pd.DataFrame(result, columns=[
                'iteration', 'seed size', 'infected', 'neighbors', 'algo', 'max propagation time', 'exp']), ignore_index=False)

        result_df = result_df.groupby(
            ['iteration', 'seed size', 'algo', 'max propagation time', 'exp']).mean().reset_index()
        result_df.to_csv(processed_filename, header=True,
                         index=False, mode='a')


final_legend_dic = {
    'nbr': 'Nbr',
    'degree': 'Degree',
    'clique': 'Clique',
    'kd': '(k, d)',
    "dk": '(d, k)'
}

processed_filename = '../output/processed_propagation_result_9a_' + dataset + '.csv'
result_df = pd.read_csv(processed_filename)
result_df['algo'] = result_df['algo'].replace(final_legend_dic)
for y_var in ['infected', 'neighbors'][:1]:

    if(dataset == "dblp"):
        result_df = result_df[result_df['seed size'] >= 300]
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.barplot(x='seed size', y=y_var, hue='algo', palette='colorblind', hue_order=[
                'Clique', 'Degree', 'Nbr', '(k, d)'], data=result_df)
    plt.xlabel('#seed from inner core', fontsize=fontsize-2)
    # plt.ylabel(y_var, fontsize=fontsize)
    plt.ylabel("Avg. spread\nper seed", fontsize=fontsize-2)
    ax.yaxis.set_major_formatter(ticker.EngFormatter())
    # ax.xaxis.set_major_formatter(ticker.EngFormatter())
    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)
    plt.grid(axis='y')
    plt.legend(loc='upper center', bbox_to_anchor=(
        0.5, 1.3), ncol=4, fancybox=False, shadow=True, fontsize=labelsize-2, columnspacing=0.8)
    # plt.title(key, fontsize=fontsize)
    plt.tight_layout()
    filename = dataset + "_" + y_var + ".pdf"
    print(filename)
    if(save):
        plt.savefig("../fig/" + filename)
        plt.show()
    plt.clf()

    # break
