#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 09:15:37 2021

@author: fra
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import sys
sys.path.append('../')

from PDE_solver import SIR_PDEroutine
from Likelihood import log_likelihood_models
"""
Usage: 

python Pde_inference.py \
    --infectious_time_distribution exponential
    --contact_interval_distribution gamma
    --lower_bounds 0.0005,0.1,1,0.1
    --upper_bounds 0.001,1,10,1
    --times_of_infections infection_times.txt
    --times_of_recovery=recovery_times.txt\
    --precision 0.1
    --output_files output
    --seed 1

"""


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--infectious_time_distribution", type = str, required = False, \
       help = "name of recovery/infectious time distribution", default = "exponential")
      
    parser.add_argument('-c','--contact_interval_distribution',  required = False, \
        help = "name of contact interval distribution",type = str, default = "exponential")

    parser.add_argument('-lb','--lower_bounds', nargs = '+', required = True, \
        help = "lower bounds for all the parameters",type = float, default = [5e-4,0.01,2,1])
            
    parser.add_argument('-ub','--upper_bounds', nargs = '+', required = True, \
        help = "upper_bounds for all the parameters",type = float, default = [1e-2,2,20,15])
    
    parser.add_argument('-ti','--infection_times',  required = False, \
        help = "path to times of infection data",type = str, default = None)
    parser.add_argument('-tr','--recovery_times',  required = False, \
        help = "path to times of recovery data",type = str, default = None)
    parser.add_argument('-p','--precision',  required = False, \
        help = "dx and dt",type = float, default = 0.2)    
    parser.add_argument('-o','--output_files',  required = False, \
        help = "name of output files",type = str, default = "result")

    parser.add_argument('-s','--seed',  required = False, \
        help = "seed for rng",type = int, default = 1)
        

    args = parser.parse_args()


    np.random.seed(args.seed)

    if args.infection_times == None and args.recovery_times==None:
        print("No data given, both infection and recovery times are null\n")
        exit(2)


    if args.contact_interval_distribution=='exponential':
        def inf_haz(u,*CIdistParms):
            beta = float(CIdistParms[0])
            return beta*np.ones_like(u)  
        def inf_distr(u,*CIdistParms):
            beta = float(CIdistParms[0])
            return stats.expon.pdf(u, 1/beta) 
        hazard_inf_par=1

    elif args.contact_interval_distribution=='gamma':
        def inf_haz(u,*CIdistParms):
            a = CIdistParms[0]
            scale = CIdistParms[1]
            tol = 1e-10
            #Basically: use de l'hopital when the ratio becomes 0/0
            #Otherwise go with definition. This regularises a lot the numerics
            x=np.zeros_like(u)
            bad_indexes=stats.gamma.cdf(u,a=a,scale=scale)>1-tol
            x[bad_indexes]=1/scale - (a-1)/u[bad_indexes]
            x[~bad_indexes]= stats.gamma.pdf(u[~bad_indexes],a=a,scale=scale)/(1- stats.gamma.cdf(u[~bad_indexes],a=a,scale=scale))
   
     
            return x   
        def inf_distr(u,*CIdistParms):
            a = CIdistParms[0]
            scale =CIdistParms[1]
            return stats.gamma.pdf(u,a=a,scale=scale)        
        hazard_inf_par=2

    elif args.contact_interval_distribution=='weibull':
        def inf_haz(u,*CIdistParms):
           shape=CIdistParms[0]
           scale=CIdistParms[1]
           return shape/scale * (u/scale)**(shape-1)  
        def inf_distr(u,*CIdistParms):
            shape = CIdistParms[0]
            scale = CIdistParms[1]
            return stats.weibull_min.pdf(u,c=shape,scale=scale)       
        hazard_inf_par=2
    
    elif args.contact_interval_distribution=='lognormal':       
        def inf_haz(u,*CIdistParms):
           s=CIdistParms[0]
           scale=CIdistParms[1]
           return stats.lognorm.pdf(u, s=s, scale=scale)/(1- stats.lognorm.cdf(u,s=s,scale=scale))
        hazard_inf_par=2
        def inf_distr(u,*CIdistParms):
            s=CIdistParms[0]
            scale=CIdistParms[1]
            return stats.lognorm.pdf(u, s=s, scale=scale)    
    else:
        print("The only distributions available are exponential, gamma, weibull,\
              lognormal. Check this code to define the functions you need \n")
        exit(1)
 
    
     #recovery/infectious period
    
    if args.infectious_time_distribution=='exponential':
        def rec_haz(u,*RecParms):
            beta = float(RecParms[0])
            return beta*np.ones_like(u)     
        def rec_distr(u,*RecParms):
            beta = float(RecParms[0])
            return stats.expon.pdf(u, 1/beta)  
        rec_parms=1
    elif args.infectious_time_distribution=='gamma':
        def rec_haz(u,*RecParms):
            a = RecParms[0]
            scale = RecParms[1]
            tol = 1e-10
            #Basically: use de l'hopital when the ratio becomes 0/0
            #Otherwise go with definition. This regularises a lot the numerics
            x=np.zeros_like(u)
            bad_indexes=stats.gamma.cdf(u,a=a,scale=scale)>1-tol
            x[bad_indexes]=1/scale - (a-1)/u[bad_indexes]
            x[~bad_indexes]= stats.gamma.pdf(u[~bad_indexes],a=a,scale=scale)/(1- stats.gamma.cdf(u[~bad_indexes],a=a,scale=scale))
   
            return x
        rec_parms=2
        
        def rec_distr(u,*RecParms):
            a = RecParms[0]
            scale = RecParms[1]
            return stats.gamma.pdf(u,a=a,scale=scale)
    elif args.infectious_time_distribution=='weibull':
        def rec_haz(u,*CIdistParms):
           shape=CIdistParms[0]
           scale=CIdistParms[1]
           return shape/scale * (u/scale)**(shape-1)  
        def rec_distr(u,*RecParms):
            shape = RecParms[0]
            scale = RecParms[1]
            return stats.weibull_min.pdf(u,c=shape,scale=scale)
        rec_parms=2

    elif args.infectious_time_distribution=='lognormal':
        
        def rec_haz(u,*RecParms):
           s=RecParms[0]
           scale=RecParms[1]
           return stats.lognorm.pdf(u, s=s, scale=scale)/(1- stats.lognorm.cdf(u,s=s,scale=scale))
        def rec_distr(u,*RecParms):
            s=RecParms[0]
            scale=RecParms[1]
            return stats.lognorm.pdf(u, s=s, scale=scale)     
        rec_parms=2

    else:
        print("The only distributions available are exponential, gamma, \
              weibull, lognormal. Check this code to define the functions you need\n")
        exit(1)    

    if args.infection_times != None:
        times_of_infections = np.loadtxt("{}".format(args.infection_times)) 
    if args.recovery_times != None:
        times_of_recovery = np.loadtxt("{}".format(args.recovery_times)) 
    
    
    
    print("Infectious time distribution: {}\nContact interval distribution: {}\n\
          times_of_infections file: {}\ntimes_of_recovery file: {}\n\
              output files prefix: {}\n".format(args.infectious_time_distribution,\
               args.contact_interval_distribution, args.infection_times, args.recovery_times, args.output_files))
    
    
    T_f = max(max(times_of_infections),max(times_of_recovery))

    grids = int(np.ceil(T_f/args.precision))
    ll=log_likelihood_models(grids,hazard_inf=inf_haz,hazard_rec=rec_haz, 
                             rec_distr = rec_distr, 
                          T=T_f, infect_times=times_of_infections,
                          recov_times=times_of_recovery,
                          hazard_inf_par=hazard_inf_par,rec_parms=rec_parms)
    result = ll.minimize_likelihood(np.array(args.lower_bounds), np.array(args.upper_bounds), maxiter=100,swarmsize=100)
         
    
    
    parameters=result[0]
    pde = SIR_PDEroutine(parameters[0], CIdist=inf_haz, CIdistParms=[parameters[1]],\
                              recovDist=rec_haz, recovDistParms=[parameters[2], parameters[3]],\
                                  nTgrid=grids, nUgrid=grids, T=T_f)
    #Initial condition in this case should be a delta in 0.
    initial_condition=np.zeros_like(pde.tgrids)
    initial_condition[0]=1
    #Solve the PDE
    S_mle,I_mle=pde.finDiffUpdate(intiial_condition=initial_condition)         
    
    np.savetxt("{}_pde.csv".format(args.output_files), np.c_[pde.tgrids,S_mle,I_mle],header="time,S,I", delimiter=',')
    np.savetxt("{}_parameters.csv".format(args.output_files), result[0],header="rho,infection_parameters,recovery_parameters", delimiter=',')
    
    x=np.linspace(0,T_f,grids)
    
    plt.figure(figsize=(6,6))
    plt.plot(pde.tgrids,I_mle, color='k')
    plt.xlabel('time')
    plt.ylabel('I(t)')
    plt.xlim(0)
    plt.ylim(0)
    plt.savefig("{}_I_fig.pdf".format(args.output_files), format='pdf')
    
    plt.close()
    
    plt.figure(figsize=(6,6))
    plt.plot(pde.tgrids,S_mle, color='k')
    plt.xlabel('time')
    plt.ylabel('S(t)')
    plt.xlim(0)
    plt.ylim(0)    
    plt.savefig("{}_S_fig.pdf".format(args.output_files), format='pdf')    
    plt.close()    
    
    t=np.linspace(0,30,600)
    plt.figure(figsize=(6,6))
    plt.plot(t,inf_distr(t,*parameters[1:1+hazard_inf_par]), color='r')
    plt.xlabel('time')
    plt.ylabel('pdf')
    plt.title("Contact interval distribution, params:{}".format(parameters[1:1+hazard_inf_par]))
    plt.xlim(0)
    plt.ylim(0)    
    plt.savefig("{}_contact_interval_fig.pdf".format(args.output_files), format='pdf')    
    plt.close()        
    plt.figure(figsize=(6,6))
    plt.plot(t,rec_distr(t,*parameters[1+hazard_inf_par:]), color='b')
    plt.xlabel('time')
    plt.ylabel('pdf')
    plt.title("Infectious period distribution, params:{}".format(parameters[1+hazard_inf_par:]))    
    plt.xlim(0)
    plt.ylim(0)    
    plt.savefig("{}_infectious_time_fig.pdf".format(args.output_files), format='pdf')    
    plt.close()        
    
    
    
    
    
    
