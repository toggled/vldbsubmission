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

# python sir_propagation_exp.py --dataset dblp --algo graph_core --num_delete -1
# python run.py -sir_exp3 1 -a graph_core -d dblp --num_delete -1
# python sir_propagation_exp.py --dataset dblp --algo naive_degree --num_delete -1
# python run.py -sir_exp3 1 -a naive_degree -d dblp --num_delete -1
# python sir_propagation_exp.py --dataset dblp --algo naive_nbr --num_delete -1
# python run.py -sir_exp3 1 -a naive_nbr -d dblp --num_delete -1






# =====================================================================
# 9(c) experiment
# =====================================================================

memlimit="3000000"
ulimit -v $memlimit

dataset="dblp"
rm ../output/propagation_result_recursive_delinner_$dataset*
for num_delete in {1000..2000..1000}
do 
	python potential_seeds.py --dataset $dataset --num_delete $num_delete
	for algo in graph_core naive_degree naive_nbr    
	do
		# python run.py -sir_exp3 1 -a $algo -d $dataset --num_delete $num_delete
		python sir_exp3.py -a $algo -d $dataset --num_delete $num_delete
	done
done


# =====================================================================
# 9(a) experiment
# =====================================================================
# rm ../output/propagation_result_9a.csv
# for i in 1000 1500 2000 3000 4000 5000 6000
# for i in 10 50 100 200 300 400 1000 1500 2000 3000
#for i in 100
#do
#    for algo in graph_core naive_degree naive_nbr    
#    do
#        python run.py --sir_9a -a $algo -d dblp --seed_size $i --max_propagation_time 100
#    done
#done
# for nd in {1000..10000..1000}
# do
# 	for i in 100
# 	do
# 	    for algo in graph_core naive_degree naive_nbr    
# 	    do
# 		python run.py --sir_9a -a $algo -d dblp --num_delete $nd --seed_size $i --max_propagation_time 100
# 	    done
# 	done
# done


# =============================================
# Plots 9(b)
# =============================================
# python plot_9b.py



# =============================================
# Statistics
# =============================================
# mkdir -p statistics
# touch statistics/neighbor.txt
# for dataset in enron dblp pref
# do 
#     python statistics_neighbor.py --dataset $dataset >> statistics/neighbor.txt
# done


# for dataset in enron dblp pref
# do 
#     python statistics_core.py --dataset $dataset
# done
