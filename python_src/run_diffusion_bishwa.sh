# python SirApplication3data.py --level 2 -a naive_nbr -d enron --num_delete 10
# python SirApplication3data.py --level 2 -a naive_degree -d enron --num_delete 10
# python SirApplication3data.py --level 2 -a graph_core -d enron --num_delete 10
# python SirApplication3data.py --level 2 -a naive_nbr -d enron --num_delete 18
# python SirApplication3data.py --level 2 -a naive_degree -d enron --num_delete 18
# python SirApplication3data.py --level 2 -a graph_core -d enron --num_delete 18


# python run.py -sir_exp3 1 -a naive_nbr -d enron --num_delete 10
# python run.py -sir_exp3 1 -a naive_degree -d enron --num_delete 10
# python run.py -sir_exp3 1 -a graph_core -d enron --num_delete 10
# python run.py -sir_exp3 1 -a naive_nbr -d enron --num_delete 18
# python run.py -sir_exp3 1 -a naive_degree -d enron --num_delete 18
# python run.py -sir_exp3 1 -a graph_core -d enron --num_delete 18


# python SirApplication3data.py --level 2 -a naive_nbr -d default --num_delete 400
# python run.py -a naive_nbr -d default -sir_exp3 1 --num_delete 400
# python SirApplication3data.py --level 2 -a naive_degree -d default --num_delete 400
# python run.py -a naive_degree -d default -sir_exp3 1 --num_delete 400
# python SirApplication3data.py --level 2 -a graph_core -d default --num_delete 400
# python run.py -a graph_core -d default -sir_exp3 1 --num_delete 400


# python SirApplication3data.py --level 2 -a naive_nbr -d default --num_delete 2
# python run.py -a naive_nbr -d default -sir_exp3 1 --num_delete 2
# python SirApplication3data.py --level 2 -a naive_degree -d default --num_delete 2
# python run.py -a naive_degree -d default -sir_exp3 1 --num_delete 2
# python SirApplication3data.py --level 2 -a graph_core -d default --num_delete 2
# python run.py -a graph_core -d default -sir_exp3 1 --num_delete 2

# python SirApplication3data.py --level 2 -a naive_nbr -d default --num_delete 4
# python run.py -a naive_nbr -d default -sir_exp3 1 --num_delete 4
# python SirApplication3data.py --level 2 -a naive_degree -d default --num_delete 4
# python run.py -a naive_degree -d default -sir_exp3 1 --num_delete 4
# python SirApplication3data.py --level 2 -a graph_core -d default --num_delete 4
# python run.py -a graph_core -d default -sir_exp3 1 --num_delete 4

# python SirApplication3data.py --level 2 -a naive_nbr -d default --num_delete 8
# python run.py -a naive_nbr -d default -sir_exp3 1 --num_delete 8
# python SirApplication3data.py --level 2 -a naive_degree -d default --num_delete 8
# python run.py -a naive_degree -d default -sir_exp3 1 --num_delete 8
# python SirApplication3data.py --level 2 -a graph_core -d default --num_delete 8
# python run.py -a graph_core -d default -sir_exp3 1 --num_delete 8



# python run.py --sir 1 --dataset dblp --max_propagation_time 1 --algo naive_nbr
# python run.py --sir 1 --dataset dblp --max_propagation_time 1 --algo naive_degree
# # python run.py --sir 1 --dataset dblp --max_propagation_time 2 --algo naive_nbr
# # python run.py --sir 1 --dataset dblp --max_propagation_time 2 --algo naive_degree

# python run.py --sir_kd 1 --dataset dblp --max_propagation_time 1 --algo naive_nbr
# python run.py --sir_kd 1 --dataset dblp --max_propagation_time 1 --algo naive_degree
# # python run.py --sir_kd 1 --dataset dblp --max_propagation_time 2 --algo naive_nbr
# # python run.py --sir_kd 1 --dataset dblp --max_propagation_time 2 --algo naive_degree


# =====================================================================
# 9(a) experiment
# =====================================================================

# rm ../output/propagation_result.csv
# time python run.py --sir 1 --dataset dblp --max_propagation_time 100 --algo naive_nbr
# time python run.py --sir 1 --dataset dblp --max_propagation_time 100 --algo naive_degree
# time python run.py --sir 1 --dataset dblp --max_propagation_time 100 --algo graph_core
# time python run.py --sir_kd 1 --dataset enron --max_propagation_time 100 --algo naive_nbr



# python SirApplication3data.py --level 2 -a naive_nbr -d enron --num_delete -1
# python run.py -sir_exp3 1 -a naive_nbr -d enron --num_delete -1
# python SirApplication3data.py --level 2 -a naive_degree -d enron --num_delete -1
# python run.py -sir_exp3 1 -a naive_degree -d enron --num_delete -1
# python SirApplication3data.py --level 2 -a graph_core -d enron --num_delete -1
# python run.py -sir_exp3 1 -a graph_core -d enron --num_delete -1



# =====================================================================
# 9(c) experiment
# =====================================================================

# python sir_propagation_exp.py --dataset dblp --algo graph_core --num_delete -1
# python run.py -sir_exp3 1 -a graph_core -d dblp --num_delete -1
# python sir_propagation_exp.py --dataset dblp --algo naive_degree --num_delete -1
# python run.py -sir_exp3 1 -a naive_degree -d dblp --num_delete -1
# python sir_propagation_exp.py --dataset dblp --algo naive_nbr --num_delete -1
# python run.py -sir_exp3 1 -a naive_nbr -d dblp --num_delete -1


# rm ../output/propagation_result_9a.csv
for algo in graph_core naive_degree naive_nbr
do
    for i in 10 50 100 200 300 400 1000 1500 2000 3000
    do
        python run.py --sir_9a -a $algo -d enron --seed_size $i --max_propagation_time 10
    done
done