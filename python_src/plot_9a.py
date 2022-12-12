import pandas as pd
import seaborn as sns
# %matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import os
from ast import literal_eval
import shutil

output_folder = '../output/'

fontsize = 18
labelsize = 14

df = pd.read_csv(output_folder + "propagation_result_9a.csv", header=None)
df.columns = ['algo', 'dataset', 'exp', 'intervention_results', 'max propagation time', 'num delete', 'p', 'result', 'seed size', 'timestep_results']
df['seed size'] = df['seed size'].astype(int)
df['max propagation time'] = df['max propagation time'].astype(int)
print(df.shape)
print(df.head(15))



import seaborn as sns
sns.set(rc={'figure.figsize': (10, 5)})
sns.set_style("ticks")
# plt.style.use('grayscale')
lw = 3
save = True

goodname_algo = {
    'naive_nbr': 'nbr',
    'naive_degree': 'degree',
    'graph_core' : 'clique'
}
df2 = df.copy()
group_list = ['dataset']
for key, item in df2.groupby(group_list, as_index=False):
    processed_filename = '../output/processed_propagation_result_9a_' + key + '_.csv'

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


    result_df = result_df.groupby(['iteration', 'seed size', 'algo', 'max propagation time' ,'exp']).mean().reset_index()
    result_df.to_csv(processed_filename, header=True, index=False, mode='a')


for key, item in df2.groupby(group_list, as_index=False):
    processed_filename = '../output/processed_propagation_result_9a_' + key + '_.csv'
    result_df = pd.read_csv(processed_filename)
    
    for y_var in ['infected', 'neighbors'][:1]:

        sns.barplot(x='seed size', y=y_var, hue='algo', data=result_df)

        plt.xlabel('#seed from inner core', fontsize=fontsize)
        # plt.ylabel(y_var, fontsize=fontsize)
        plt.ylabel("Avg. spread per seed", fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(loc='best', fontsize=fontsize)
        plt.title(key, fontsize=fontsize)
        plt.tight_layout()
        filename = key + "_" + y_var + ".pdf"
        print(filename)
        if(save):
            plt.savefig("../fig/" + filename)
            # plt.show()
        plt.clf()

    # break
