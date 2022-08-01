# python -W ignore -u run.py --algo naive_nbr --dataset default --sir
# python -W ignore -u run.py --algo naive_nbr --dataset default --sir_exp2
# python -W ignore -u run.py --algo naive_nbr --dataset default --sir_exp3
# python SirApplication3data.py -d default -a naive_nbr -l 2
# python -W ignore -u run.py --algo naive_nbr --dataset default --sir_exp3_explanation
# python -W ignore -u run.py --algo naive_nbr --dataset default --sir_exp3_explanation_splen

python -W ignore -u run.py --algo naive_nbr --dataset enron --sir
python -W ignore -u run.py --algo naive_nbr --dataset enron --sir_exp2
python -W ignore -u run.py --algo naive_nbr --dataset enron --sir_exp3
python SirApplication3data.py -d enron -a naive_nbr -l 5
python -W ignore -u run.py --algo naive_nbr --dataset enron --sir_exp3_explanation
python -W ignore -u run.py --algo naive_nbr --dataset enron --sir_exp3_explanation_splen
