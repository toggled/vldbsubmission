import os
import argparse
import pickle
parser = argparse.ArgumentParser()
parser.add_argument("--s", help="token of dir", default="Local-core(P)_enron", type=str)
parser.add_argument("--t", help="token of dir", default="E-Peel_enron", type=str)
# parser.add_argument("--all", action='store_true')

args = parser.parse_args()
def read(filename):
    print("reading file: ",filename)
    _core = {}
    with open(filename,'r') as f:
        for line in f.readlines():
            x,y = line.split(',')
            x,y = x.strip(),y.strip()
            _core[x] = int(y)
    return _core 

def singletest(source_fname,target_fname):
    core_compared,core_base = read(source_fname),read(target_fname)
    # print(core_compared)
    # print(core_base)
    assert len(core_base) == len(core_compared), "Two returned cores do not have same length: " + str(len(core_base)) + " != " + str(len(core_compared))
    for v in core_base:
        assert v in core_compared, str(v) + " is not in core_compared"
        assert core_base[v] == core_compared[v], str(v)+" :Output core is different in " + str(core_base[v]) + " & " + str(core_compared[v])

    # core_compared contains in core_base
    for v in core_compared:
        assert v in core_base, str(v) + " is not in core_base"
        assert core_base[v] == core_compared[v], str(v)+" :Output core is different in " + str(core_base[v]) + " & " + str(core_compared[v])

    print("\n Passed")

singletest(args.s, args.t)
# source_fname = "output/core_"+args.s+".csv"
# target_fname = "output/core_"+args.t+".csv"

# singletest(source_fname,target_fname)

# else:
#     datasets = ["enron"]
#     # datasets = ['bin_1']
#     # datasets = ["syn","enron" ,"bin_1", "bin_2", "bin_4", "bin_5"]
#     # datasets = ["syn","enron" ,"bin_1", "bin_2", "bin_4", "bin_5", "contact", "congress"]
#     # datasets = ["syn","enron" ,"bin_1", "bin_2", "bin_4", "bin_5"]
#     # algos = ['Peel']
#     # algos = ['Peel','E-Peel']
#     algos = ['Local-core']
#     # algos = ['Local-core(P2)']
#     # algos =['Local-core-OPTIII(P2)']
#     # algos = ["Local-core","Local-core-OPTI", "Local-core-OPTII" ,"Local-core-OPTIII" ,"Local-core-OPTIV"]

#     for data in datasets:
#     #     target = "E-Peel_"+data
#     #     for alg in algos:
#     #         source_fname = "output/core_"+alg+"_"+data+".csv"
#     #         target_fname = "output/core_"+target+".csv"
#     #         singletest(source_fname,target_fname)

#         pyfname = "/Users/nus/hg-core-decomposition/tests/tmp/" + data + ".pkl"
#         with open(pyfname, 'rb') as handle:
#             hgDecompose = pickle.load(handle)
#         core_base = hgDecompose.core
#         # print(core_base)
#         for alg in algos:
#             source_fname = "output/core_"+alg+"_"+data+".csv"
#             core_compared = read(source_fname)
#             # print(core_compared)
#             try:
#                 assert len(core_base) == len(core_compared), "Two returned cores do not have same length: " + str(len(core_base)) + " != " + str(len(core_compared))
#                 for v in core_base:
#                     assert v in core_compared, str(v) + " is not in core_compared"
#                     assert core_base[v] == core_compared[v], str(v)+" :Output core is different in " + str(core_base[v]) + " & " + str(core_compared[v])

#                 # core_compared contains in core_base
#                 for v in core_compared:
#                     assert v in core_base, str(v) + " is not in core_base"
#                     assert core_base[v] == core_compared[v], str(v)+" :Output core is different in " + str(core_base[v]) + " & " + str(core_compared[v])
#                 print("\n Passed")
#             except Exception as e:
#                 print('\n failed')
#                 print(e)
            
