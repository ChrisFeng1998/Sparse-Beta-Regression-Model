# A Sparse Beta Regression Model for Network Analysis

Code for reproducing the results shown in the paper 'A Sparse Beta Regression Model for Network Analysis'

## Simulation results

The code for reproducing the simulation results are stored in the subdirectories `model_i, i=1,2,3,4,5,6`.

Each of those subdirectories contains the following:

- `run_simulation.R` contains the code for running the simulation, see instructions below.
- the `simulation_results` subdirectory contains the raw output from the simulations
- `generate_summary.Rmd` contains to code for producing the simulations summaries presented in the paper, such as confidence intervals and error plots.


### Re-running the simulations

The simulation data in the subdirectories `model_i/simulation_results` can be reproduced by running the `run_simulation.R` files. The file will look for the environment variable `SGE_TASK_ID` and uses its value for setting the seed. Running the file once will do 100 simulation runs each for each sample size. To produce the full 1000 simulation runs, the file is meant to be run 10 times, each time with a different seed.

## Lawyer network

The results from analyzing Lazega's lawyer data can be reproduced by running `lawyer_friendship_network.R`. Before running the code, please download the necessary data from the [original data owner's website](https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm) and place it in the `LazegaLawyers` subdirectory.

## World trade netwrok

The results from the analysis of the world trade network can be reproduced by running `trade_network.R`. Before running the code, please download the necessary data from the [original data owner's website](http://personal.lse.ac.uk/tenreyro/LGW.html) and place it in the `trade_network` subdirectory
