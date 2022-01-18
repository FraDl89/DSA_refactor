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
└── Synthetic_data
    ├── Code
    ├── Data
    └── synthetic_Gamma

In the main folder, you can find the necessary files to run any other code present in this repo. in the file requirements.txt, all the libraries need to reproduce the results are listed (pip install requirements.txt should be sufficient to get everything in place).

To start working with the code, it is recommended to clone the repo locally and run the codes contained in Examples.

Examples/pde_example.py produces the solution of a pde. To use it, you can type in a terminal
	python pde_example.py -i gamma -c exponential -r 0.002 -t 130 -d 0.1 -ip 0.666 13.5 -cp 0.26 -o Results/pde_example

This solves the pde assuming that the infectious time distribution (-i) is gamma with parameters (-ip) shape=0.666, scale=13.5, the contact interval distribution (-c) is exponential with scale parameter (-cp) 0.26, and the initial condition (-r) is rho=0.002. The mesh (-d) both in dx and dt is of size 0.1, and the PDE is solved until the final time (-t) 130. Finally, all the outputs produced are saved in (-o) the local folder "Results/" with prefix "result". Make sure to have the Results folder inside Examples. Distributions available at the moment are 'exponential', 'gamma', 'lognormal', 'weibull'. So, for instance, you can run the same algorithm with python pde_example.py -i lognormal -c gamma -r0.02 -t 130 -d 0.1 -ip 1 1 -cp 1 1 -o result. The order of the parameters can be found by looking at the scipy.stats documentation for the relative distribution.

Examples/python Pde_inference.py lets you run inference on parameters. This uses two files (recovery_times.txt, infection_times.txt) that consist of one array each containing the times of recovery/infections. At least one file is needed. To run this code, you can type in a terminal

	python Pde_inference.py --infectious_time_distribution gamma --contact_interval_distribution exponential --lower_bounds 0.0005 0.1 1 0.1 --upper_bounds 0.005 1 10 1 --infection_times infection_times.txt --recovery_times recovery_times.txt --seed 1 --output Results/inference

This tells the likelihood scheme to infer a gamma distribution as infectious time distribution, an exponential distribution as contact interval distribution, in the region defined by lower and upper bounds for (rho,scale[exponential],shape[gamma],scale[gamma]) as indicated by --lower_bounds and --upper_bounds, respectively. Since the solver relies on a stochastic particle swarm algorithm, a seed is needed for reproducibility. All the results, which consist of plots of the solution, plots of the infectious period and contact interval distributions, and the numerical value of the solution, are saved in the local folder Results, with suffix "output". Distributions available at the moment are 'exponential', 'gamma', 'lognormal', 'weibull'. So, for instance, you can run the same algorithm with 
python Pde_inference.py --infectious_time_distribution gamma --contact_interval_distribution weibull --lower_bounds 0.0005 0.1, 0.1 1 0.    1 --upper_bounds 0.005 1 1 10 1 --infection_times infection_times.txt --recovery_times recovery_times.txt --seed 1 --output Results/weibull_gamma

Note, for distributions depending on more than one parameter, the order for the lower and upper bounds is always as follows: (shape,scale). In general, I tried to follow scipy.stats documentation for the relative distributions.


Examples/epidemics.py is a toy code that you can read through to understand how to generate some data, run the likelihood inference, and solve the pde.
You can copy the relevant bits of the code in a new file if you want to play with it. Since most of the code relies on defining hazard functions and probability distributions, if you want to run the PDE with specific distributions, you have to code them in Python.

The rest of the repo is as follows. In the main folder you can find the functions that are at the core of this repo.
	1) PDE_solver.py, which contains the pde solver, along with a fully worked out example
	2) Likelihood.py, which encodes all the likelihoods as defined in the paper, and it is called when maximising likelihoods in the numerics examples
	3) Selkealgo.py is a version of the Sellke algorithm, used to generate epidemic curves from non-markovian epidemic models with the Sellke construction, along with a fully worked out example, along with a fully worked out example.

The code to reproduce numerics can be found under the folders named Synthetic_data, India, and foot_and_mouth. Results and Data are already stored under the subfolders named Data, except for India data which generates files too big for github to handle. 

Figures in the paper can be reproduced by running the codes inside the Figures folder in R.


If you use or take inspiration from the code of this repo, please cite the paper "Dynamic Survival Analysis for non-Markovian Epidemic Models", by Di Lauro, KhudaBukhsh, Kiss, Kenah, Rempala. See also the Licence file.
