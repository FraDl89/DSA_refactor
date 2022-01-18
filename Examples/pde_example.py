#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright 2020-2022 Francesco Di Lauro. All Rights Reserved.
See Licence file for details.
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
    --infectious_time_params 0.6666 13.5
    --rho 0.02
    --contact_interval_params 0.2
    --dx 0.1
    --T_f 120
    --output_files output
"""

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--infectious_time_distribution", type = str, required = False, \
       help = "name of recovery/infectious time distribution", default = "exponential")
      
    parser.add_argument('-c','--contact_interval_distribution',  required = False, \
        help = "name of contact interval distribution",type = str, default = "exponential")

    parser.add_argument('-r','--rho',  required = True, \
        help = "initial condition",type = float)
    
    parser.add_argument('-t','--T_f',  required = True, \
        help = "final time",type = float)
            
    parser.add_argument('-ip','--infectious_time_params', nargs = '+', required = True, \
        help = "infectious time distribution parameters",type = float)
            
    parser.add_argument('-cp','--contact_interval_params', nargs = '+', required = True, \
        help = "contact interval distribution parameters",type = float)        
    
    parser.add_argument('-d','--dx',  required = False, \
        help = "dx and dt",type = float, default = 0.2)            
    parser.add_argument('-o','--output_files',  required = False, \
        help = "name of output files",type = str, default = "result")
       
    args = parser.parse_args()
        
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
            
            x = np.where(stats.gamma.cdf(u,a=a,scale=scale)>1-tol,
                         1/scale - (a-1)/u,
                         stats.gamma.pdf(u,a=a,scale=scale)/(1- stats.gamma.cdf(u,a=a,scale=scale)))
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
            
            x = np.where(stats.gamma.cdf(u,a=a,scale=scale)>1-tol,
                         1/scale - (a-1)/u,
                         stats.gamma.pdf(u,a=a,scale=scale)/(1- stats.gamma.cdf(u,a=a,scale=scale)))
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
        
    grids = int(args.T_f/args.dx)
    pde = SIR_PDEroutine(args.rho, CIdist=inf_haz, CIdistParms=args.contact_interval_params,\
                              recovDist=rec_haz, recovDistParms=args.infectious_time_params,\
                                  nTgrid=grids, nUgrid=grids, T=args.T_f)    
    
    initial_condition=np.zeros_like(pde.tgrids)
    initial_condition[0]=1
    #Solve the PDE
    S_mle,I_mle=pde.finDiffUpdate(intiial_condition=initial_condition)  
    
    np.savetxt("pde_{}.csv".format(args.output_files), np.c_[pde.tgrids,S_mle,I_mle],header="time,S,I", delimiter=',')

    
    plt.figure(figsize=(6,6))
    plt.plot(pde.tgrids,I_mle, color='k')
    plt.xlabel('time')
    plt.ylabel('I(t)')
    plt.xlim(0)
    plt.ylim(0)
    plt.savefig("I_fig_{}.pdf".format(args.output_files), format='pdf')
    
    plt.close()
    
    plt.figure(figsize=(6,6))
    plt.plot(pde.tgrids,S_mle, color='k')
    plt.xlabel('time')
    plt.ylabel('S(t)')
    plt.xlim(0)
    plt.ylim(0)    
    plt.savefig("S_fig_{}.pdf".format(args.output_files), format='pdf')    
    plt.close()    
    
        
        
        
        
        
