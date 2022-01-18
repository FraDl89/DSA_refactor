#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2020-2022 Francesco Di Lauro. All Rights Reserved.
See Licence file for details.
"""
import numpy as np
import sys
import scipy.stats as stats
from scipy.optimize import minimize

from pyDOE import lhs
sys.path.append('../../')

from Likelihood import log_likelihood_models


if __name__ == "__main__":
    '''
    Test on synthetic data to evaluate the performance of l_R2
    '''    
    N=10000
    T_f=80
    seed = int(sys.argv[1])   #Seed should be the name of the key in the dataset,
    #seed = 0

    file = np.load('Outputs/synthetic_Gamma/Gamma_distr_seed%d.npy'%seed)
    

    np.random.seed(seed)
    time = file[0]
    I = file[1]
    S = file[2]
    
    index = np.where(np.diff(I)<0)
    time_likelihood = time[index[0]+1]
    recovery_times = time_likelihood[time_likelihood<T_f]
    
    indextosample = np.where(N-S>0)
    time_for_likelihood =time[indextosample]
    I_tosample = I[indextosample]    
    
    index_I = np.diff(I_tosample,prepend=[I[0]]) #prepend I[0] so that you don't fit to the first point and you get a vector of the correct length
            
    infection_times = time_for_likelihood[np.where(index_I>0)]   #These are the infection times
    

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


    ll=log_likelihood_models(10000,hazard_inf=inf_distr,hazard_rec=rec_haz, 
                             rec_distr = rec_distr, 
                          T=T_f, infect_times=infection_times,recov_times=recovery_times,hazard_inf_par=1,rec_parms=2)
    result = ll.minimize_likelihood(np.array([5e-4,0.01,2,1]), np.array([1e-2,2,20,15]))
    
    print(result[0])

    
        
    
    
    
    
    
    
