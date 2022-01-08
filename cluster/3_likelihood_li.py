#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 13:48:27 2020

@author: francesco
"""
import numpy as np
import sys
from Likelihood import log_likelihood_models
import scipy.stats as stats
from PDE_solver import SIR_PDEroutine
from scipy.optimize import minimize

from pyDOE import lhs
import gc


if __name__ == "__main__":
    '''
    Test on synthetic data to evaluate the performance of l_I i.e. only times of infection
    '''
    N=10000    
    T_f=80
    seed = int(sys.argv[1])   #Seed should be the name of the key in the dataset,
    #seed = 0

    file = np.load('Outputs/synthetic_Gamma/Gamma_distr_seed%d.npy'%seed)
    results = []
    

    np.random.seed(seed)
    time = file[0]
    I = file[1]
    S = file[2]
    
    
    indextosample = np.where(N-S>0)
    time_for_likelihood =time[indextosample]
    I_tosample = I[indextosample]    
    
    index_I = np.diff(I_tosample,prepend=[I[0]]) #prepend I[0] so that you don't fit to the first point and you get a vector of the correct length
    
    
    time_to_fit = time_for_likelihood[np.where(index_I>0)]   #These are the infection times
    #Modify according to which distribution you are using!
    #Remember, you need hazard functions, not distributions
        
    def rec_haz(u, *recovDistParams):    
        a = float(recovDistParams[0])**2/float(recovDistParams[1])
        scale = float(recovDistParams[1])/float(recovDistParams[0])  
        tol = 1e-10
        #Basically: use de l'hopital when the ratio becomes 0/0
        #Otherwise go with definition. This regularises a lot the numerics
        
        x = np.where(stats.gamma.cdf(u,a=a,scale=scale)>1-tol,
                     1/scale - (a-1)/u,
                     stats.gamma.pdf(u,a=a,scale=scale)/(1- stats.gamma.cdf(u,a=a,scale=scale)))
        return x
        
    def rec_distr(u, *recovDistParams):
        a = float(recovDistParams[0])**2/float(recovDistParams[1])
        scale = float(recovDistParams[1])/float(recovDistParams[0])   
       
        return stats.gamma.pdf(u,a=a,scale=scale)
   
    def inf_distr(u,*CIdistParms):
        beta = float(CIdistParms[0])
        return beta*np.ones_like(u)          
    
    import time as tempo
    start = tempo.time()    
    ll=log_likelihood_models(10000,hazard_inf=inf_distr,hazard_rec=rec_haz, rec_distr = rec_distr, 
                              T=T_f, infect_times=time_to_fit,hazard_inf_par=1,rec_parms=2)

    result = ll.minimize_likelihood(np.array([5e-4,0.01,2,1]), np.array([1e-2,2,20,15]))
    end = tempo.time()
    print(*result.x)
    #del fun
    del result
    del ll
    gc.collect()


    
    
    
    
