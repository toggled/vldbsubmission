# Neighborhood-based Hypergraph core decomposition
This repository is the implementation of a PVLDB 2023 paper: "Neighborhood-based Hypergraph core decomposition". 

![](intro.PNG)

Please cite our <b>[extended arXiv version](https://arxiv.org/abs/2301.06426)</b> as:
```
@article{arafat2023neighborhoodbased,
      title={Neighborhood-based Hypergraph Core Decomposition}, 
      author={Naheed Anjum Arafat and Arijit Khan and Arpit Kumar Rai and Bishwamittra Ghosh},
      year={2023},
      eprint={2301.06426},
      archivePrefix={arXiv},
      primaryClass={cs.SI},
      url={https://arxiv.org/abs/2301.06426}
}
```
## Datasets
Our datasets can be downloaded from https://drive.google.com/file/d/12cvz-XtfQUbmj-gqlT9z4DKGMuhqLGjs/view?usp=sharing
- data/dataset/real : real-world datasets (UPDATE 6/10/2023: The earlier version of Aminer dataset had repeated nodes in some hyperedges. A clean version of Aminer has been uploaded) 
- data/dataset/synthetic : synthetically generated datasets
- data/dataset/protein: CORUM Protein complex hypergraph
- data/dataset/meetup: Nashville meetup dataset 
<!-- - data/dataset/kenneth_lay: Ego hypergraph of Kenneth lay (Enron's founder & CEO) -->

### Format of a hypergraph file (.hyp)
Each line is a hyperedge consisting of comma separated node ids. 

## Notebooks for plots and case-studies
- notebooks/:
  - PlotGeneration.ipynb: Notebook for plots in the paper
  - python_src/prev_9a.ipynb: Notebook for Diffusion Application plots.
  - densest_subhypergraph.ipynb: Effectiveness analysis & plots for volume densests and degree densest subhypergraphs.
  - casestudyI.ipynb: Notebook for case study I (CORUM Protein complex)
  - meetupCase.ipynb: Notebook for case study II (Meetup hypergraph)

## source codes
- cpp_src/: 
  - C++ implementation of the E-Peel, Peel, Local-core, Local-core with the optimizations, (k,d)-core, 3 baselines (clique-graph, dist-2 bipartite graph & degree-based algorithms) (`algorithms.cpp`)
- casestudycpp/:
  - C++ implementation of Densest subhypergraph extraction (`densest_subhypergraph.cpp`)
- par_src/: 
  - OpenMP implementation of Parallel local-core (`parallel_localcoreinitOpt2.cpp`)
- python_src/: Python code to run Application I (Diffusion).

## Output
output/: Output folder that contains algorithm outputs.

## How to run:
- Sequential Algorithms: Peel, E-Peel, Local-core, Local-core+OptI, Local-core+OptI+II,Local-core+OptI+II+III and Local-core(opt): 
  - `cd cpp_src` 
  - `bash run.sh` // Deactivate(Activate) logging with `log=0`(`log=1`), set #runs of an alg. with var `it`
- Parallel Local-core algorithm:
  - `cd par_src`
  - `bash hgrun.sh` // Deactivate(Activate) load-balancing by setting lb=0(lb=1) 
- Application I (Diffusion):
  - `cd python_src`
  - `bash run_diffusion.sh`

## C++ compiler requirements:
- C++11
- OpenMP 3.5 

## Python and notebook requirements:
- python 3.8 or above
- dependencies: 
  - disjoint-set (`pip install disjoint-set`)
  - hypernetx (`pip install hypernetx`)
  - pandas (`pip install pandas`)
  - numpy (`pip install numpy`)
  - seaborn (`pip install seaborn`)
  - matplotlib (`pip install matplotlib`)

<!-- ## Reproducing the plots:
  - Install python dependencies
  - Download `output.zip` file from https://drive.google.com/file/d/1-1oZc-ajChoZkRu9C-4Vbw6VFGqB-YtO/view?usp=sharing and extract in the repository folder.
  - Run the notebooks plots.ipynb, python_src/propationPlots.ipynb, densest_subhypergraph.ipynb, casestudyI.ipynb and casestudyII.ipynb -->
