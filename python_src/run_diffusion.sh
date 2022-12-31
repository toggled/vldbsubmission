# =====================================================================
# 9(c) experiment
# =====================================================================

#memlimit="3000000"
#ulimit -v $memlimit

dataset="dblp"
rm ../output/propagation_result_recursive_delinner_$dataset*
for num_delete in {1000..7000..1000}
do 
	python potential_seeds.py --dataset $dataset --num_delete $num_delete
 	for algo in graph_core naive_degree naive_nbr    
 	do
 		python sir_exp3.py -a $algo -d $dataset --num_delete $num_delete --seed_size 500
 	done
done


# =====================================================================
# 9(a) experiment
# =====================================================================
rm ../output/propagation_result_9a.csv

dataset="dblp"
for i in 300 400 1000 1500 2000 3000
do
   for algo in graph_core naive_degree naive_nbr    
   do
       python sir_exp.py --sir_9a -a $algo -d $dataset --seed_size $i --max_propagation_time 100
   done
   
   # For kd core
   python sir_exp.py --sir_kd -d $dataset --seed_size $i --max_propagation_time 100

   # For dk core
   #python sir_exp.py --sir_dk -d $dataset --seed_size $i --max_propagation_time 100

done
python plot_9a.py

# =============================================
# Plots 9(b)
# =============================================
python plot_9b.py



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
