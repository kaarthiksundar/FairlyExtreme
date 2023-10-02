# Source code for the paper "Fairly Extreme: Minimizing Outages Equitably" 

This repository contains the source code for all the case studies in for the manuscript "Fairly Extreme: Minimizing Outages Equitably" by Kaarthik Sundar, Deepjyoti Deka and Russell Bent. 

The Project.toml files contains all the Julia packages required for running the code. The two solvers used in the runs include CPLEX and Ipopt. CPLEX is a commercial MILP/MISOCP solver and hence, users require a license for the same. 

The *.sh files contain the scripts for all the case studies in the paper and the output/dc-opf-ls.zip file contains all the solutions of the runs. Also, the code to generate all the plots in the paper are written in Python and can be found in output/plots as a single notebook file. It uses all the results in the zip file to generate the plots. 

Any help in re-running or modifying the code can be directed to Kaarthik Sundar through an issue in this repository or by e-mail. 