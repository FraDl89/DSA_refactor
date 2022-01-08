This repo contains all the code to reproduce the results from the paper ...

PDE_solver.py contains the pde solver, along with a fully worked out example
Likelihood.py encodes all the likelihoods as defined in the paper, and it is called when maximising likelihoods in the numerics examples
Selkealgo.py is a version of the Sellke algorithm, used to generate epidemic curves from non-markovian epidemic models with the Sellke construction

The code to reproduce numerics can be found under the folders Syntetic_data, India, foot_and_mouth. Results and Data are already stored under the subfolders named Data, except for India data which generates files too big for github to handle. 

Figures in the paper can be reproduced by running the codes inside the Figures folder in R.

Finally, cluster folder includes some instruction to run synthetic data inference on a cluster. It uses outdated versions of the Likelihood and PDE solver, but it should work fine.

GNU General Public License?
