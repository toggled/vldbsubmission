# This repository contains the code and data used in the experiments of our paper titled: **Neighborhood-based Hypergraph core decomposition**
## Datasets
Our datasets can be downloaded from https://drive.google.com/file/d/1--B6vLUL1DF4BI6UH3d0UyXyZQJOu6Nb/view?usp=sharing 
- data/dataset/real : real-world datasets
- data/dataset/synthetic : synthetically generated datasets
- data/dataset/protein: CORUM Protein complex hypergraph
- data/dataset/kenneth_lay: Ego hypergraph of Kenneth lay (Enron's founder & CEO)

### Format of a hypergraph file (.hyp)
Each line is a hyperedge consisting of comma separated node ids. 

## Notebooks for plots and case-studies
- plots.ipynb: Notebook for plots in the paper
- densest_subhypergraph.ipynb: Effectiveness analysis & plots for volume densests and degree densest subhypergraphs.
- casestudyI.ipynb: Notebook for case study I (CORUM Protein complex)
- casestudyII.ipynb: Notebook for case study II (Kenneth lay ego hypergraph)

## source codes
- cpp_src/: C++ implementation of the E-Peel, Peel, Local-core, Local-core with the optimizations.
- par_src/: Openmp implementation of Parallel local-core.
- python_src/: Python code to run Application I (Diffusion).
- densest_src/: C++ implementation of densest subhypergraph extraction.

## Output
output/ : Output folder that contains algorithm outputs.

## How to run:
- Sequential Algorithms: Peel, E-Peel, Local-core, Local-core+OptI, Local-core+OptI+II,Local-core+OptI+II+III and Local-core(opt): 
  - `cd cpp_src` 
  - `bash run.sh` // Deactivate(Activate) logging with `log=0`(`log=1`), set #runs of an alg. with var `it`
- Parallel algorithm: Parallel implementation of Local-core:
  - `cd par_src`
  - `bash hgrun.sh` // Deactivate(Activate) load-balancing by setting lb=0(lb=1) 
- Application I (Diffusion):
  - TODO
- Application II (Densest subhypergraph extraction)
  - `cd cpp_src`
  - `bash denseSrun.sh`
