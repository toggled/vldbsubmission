#!/bin/bash
python run.py -a naive_nbr -d enron -sir 1
python3 run.py -a naive_nbr -d enron -sir_exp2 1
python3 run.py -a naive_nbr -d enron -sir_exp3 1
python3 run.py -a naive_degree -d enron -sir_exp3 1
python3 run.py -a graph_core -d enron -sir_exp3 1
python3 SirApplication3data.py -d enron -a naive_nbr -l 5
python3 run.py -a naive_nbr -d enron -sir_exp3_explanation 1
python3 run.py -a naive_degree -d enron -sir_exp3_explanation 1
python3 run.py -a graph_core -d enron -sir_exp3_explanation 1
python3 run.py -a naive_nbr -d enron -sir_exp3_explanation_splen 1
python3 run.py -a naive_degree -d enron -sir_exp3_explanation_splen 1
python3 run.py -a graph_core -d enron -sir_exp3_explanation_splen 1