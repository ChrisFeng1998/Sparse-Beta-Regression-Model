# A Sparse Beta Regression Model for Network Analysis

Code for reproducing the results shown in the paper 'A Sparse Beta Regression Model for Network Analysis'

## Simulation results

### Sparse beta regression model
The code for reproducing the simulation results for the sparse beta regression model is stored in the directories `model_i, i=1,2,3,4,5,6`.

Each of those directories contains the following:

- `run_simulation.R` contains the code for running the simulation, see instructions below.  The file will look for the environment variable `SLURM_JOB_ID` and use its value for setting the seed. Running the file once will do 100 simulation runs for each sample size. To produce the full 1000 simulation runs, the file is meant to be run 10 times, each time with a different seed.
- the `simulation_results` subdirectory contains the raw output from the simulations
- `generate_summary.Rmd` contains code for producing the simulation summaries presented in the paper, such as confidence intervals and error plots.


### ERC model

The code for reproducing the simulation results for the ERC model is stored in the subdirectories `erc/erc_i, i=1,2,3`.

Each of those subdirectories contains the following:

- `ER-C.R` contains the code for running the simulation, see instructions below. As with the other simulation results, the script will look for the `SLURM_JOB_ID` environment variable for setting the seed. Running the script once will run 20 simulations for each sample size and to reproduce the full simulation results the script is meant to be run 50 times, each time with a different seed.
- `ER-C_generate_summary.Rmd` contains code for producing the simulation summaries presented in the paper. Knitting the file will produce an `ER-C_generate_summary.pdf` and `ER-C_generate_summary.tex` file, containing the table with the values for the confidence intervals.

## Data Analysis

### Lawyer network

The results from analyzing Lazega's lawyer data can be reproduced by running `lawyer_friendship_network.R`. Before running the code, please download the necessary data from the [original data owner's website](https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm).

### World trade network

The results from the analysis of the world trade network can be reproduced by running `trade_network.R`. Before running the code, please download the necessary data from the [original data owner's website](http://personal.lse.ac.uk/tenreyro/LGW.html).
