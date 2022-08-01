
import os
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("--thread", help="index of thread", default=-1, type=int)
parser.add_argument("--max_thread", help="maximum number of thread", default=1, type=int)
args = parser.parse_args()

# algo_list = ['naive_nbr', 'naive_degree']
algo_list = ['naive_nbr']
# dataset_list = ['bin_1', 'bin_2', 'bin_4', 'bin_5', 'enron',  'contact', 'congress']
# dataset_list = ['congress']
dataset_list = ['enron']
# exps = ['sir', 'sir_exp2', 'sir_exp3']
exps = ['sir_exp3']

# all combination of experiments
configurations = []
for dataset in dataset_list:
    for algo in algo_list:
        for exp in exps:
            configurations.append((algo, dataset, exp))


# print(len(configurations))
# distributing among threads
for i, configuration in enumerate(configurations):
    algo, dataset, exp = configuration
    if(i%args.max_thread == args.thread or args.thread == -1):
        cmd = "python -W ignore -u run.py" + \
              " --algo " + algo + \
              " --dataset " + dataset + \
              " --" + exp
        print(cmd) 
        # os.system(cmd) 

# TO DO: ignore assertion -O
os.system('python SirApplication3data.py -d default -a naive_nbr -l 5')