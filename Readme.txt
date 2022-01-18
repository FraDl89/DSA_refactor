This repo contains all the code to reproduce the results from the paper "Dynamic Survival Analysis for non-Markovian Epidemic Models".

This is the structure of the Repo.
.
├── Examples
├── Figures
│   ├── foot_and_mouth_conf_int
│   ├── India
│   ├── __pycache__
│   └── Synthetic
├── foot_and_mouth
│   ├── Bootstrap
│   ├── Code
│   └── Data
├── India
│   ├── Bootstrap
│   ├── Code
│   └── Data
├── cluster
└── Synthetic_data
    ├── Code
    ├── Data
    └── synthetic_Gamma

To start working with the code, it is recommended to clone the repo locally and run the codes contained in Examples.
Examples/pde_example.py produces the solution of a pde. To use it, you can type in a terminal
	python pde_example.py -i gamma -c exponential -r 0.002 -t 130 -d 0.1 -ip 0.666 13.5 -cp 0.26 -o result
This solves the pde assuming that the infectious time distribution is gamma with parameters shape=0.666, scale=13.5, the contact interval distribution is exponential with scale parameter 0.26, and the initial condition is rho=0.002. The grid (both in dx and dt) is of length 0.1, and the PDE is solved until T_f=130. Finally, all the outputs produced are saved with suffix "result". Distributions available at the moment are 'exponential', 'gamma', 'lognormal', 'weibull'. So, for instance, you can run the same algorithm with python pde_example.py -i lognormal -c gamma -r0.02 -t 130 -d 0.1 -ip 1 1 -cp 1 1 -o result. The order of the parameters can be found by looking at the scipy.stats documentation for the relative distribution.

Examples/python Pde_inference.py lets you run inference on parameters. This uses two files (recovery_times.txt, infection_times.txt) that consist of one array each containing the times of recovery/infections. At least one file is needed. To use it, you can type in a termina
	python Pde_inference.py --infectious_time_distribution gamma --contact_interval_distribution exponential --lower_bounds 0.0005 0.1 1 0.1 --upper_bounds 0.005 1 10 1 --infection_times infection_times.txt --recovery_times recovery_times.txt --seed 1 --output output
This tells the likelihood scheme to infer a gamma distribution as infectious time distribution, an exponential distribution as contact interval distribution, in the region defined by lower and upper bounds for (rho,scale[exponential],shape[gamma],scale[gamma]) as indicated by --lower_bounds and --upper_bounds, respectively. Since the solver relies on a stochastic particle swarm algorithm, a seed is needed for reproducibility. All the results, which consist of plots of the solution, plots of the infectious period and contact interval distributions, and the numerical value of the solution, are saved with the suffix "output". Distributions available at the moment are 'exp    onential', 'gamma', 'lognormal', 'weibull'. So, for instance, you can run the same algorithm with python pde_example.py -i lognormal -c gamma -    r0.02 -t 130 -d 0.1 -ip 1 1 -cp 1 1 -o result. The order of the parameters can be found by looking at the scipy.stats documentation for the relative distribution.


Examples/epidemics.py is a toy code that you can read through to understand how to generate some data, run the likelihood inference, and solve the pde.
You can copy the relevant bits of the code in a new file if you want to play with it. Since most of the code relies on defining hazard functions and probability distributions, if you want to run the PDE with specific distributions, you have to code them in Python.

The rest of the repo is as follows. In the main folder you can find
	1) PDE_solver.py, which contains the pde solver, along with a fully worked out example
	2) Likelihood.py, which encodes all the likelihoods as defined in the paper, and it is called when maximising likelihoods in the numerics examples
	3) Selkealgo.py is a version of the Sellke algorithm, used to generate epidemic curves from non-markovian epidemic models with the Sellke construction

The code to reproduce numerics can be found under the folders Syntetic_data, India, foot_and_mouth. Results and Data are already stored under the subfolders named Data, except for India data which generates files too big for github to handle. 

Figures in the paper can be reproduced by running the codes inside the Figures folder in R.


This repo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

If you use or take inspiration from the code of this repo, please cite the paper ""
