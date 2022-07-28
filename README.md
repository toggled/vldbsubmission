# This repository contains the code and data used in the experiments of our paper titled:
"Neighborhood-based Hypergraph core decomposition"
## Data
Our datasets can be downloaded from https://drive.google.com/file/d/1--B6vLUL1DF4BI6UH3d0UyXyZQJOu6Nb/view?usp=sharing 
- data/dataset/real : real-world datasets
- data/dataset/synthetic : synthetically generated datasets
- data/dataset/protein: CORUM Protein complex dataset
- data/dataset/kenneth_lay: Ego hypergraph of Kenneth lay (enron)

## Notebooks
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
- Peel, E-Peel, Local-core and Local-core(opt): `cd cpp_src` `bash cpp_src/run.sh`
