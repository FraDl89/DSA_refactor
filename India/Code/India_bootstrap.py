#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 13:35:06 2021

@author: fra
"""
import numpy as np
import scipy.stats as stats
import random
import sys


sys.path.append('../../')
from Likelihood import log_likelihood_models
from PDE_solver import SIR_PDEroutine


if __name__=="__main__":
    
    seed =  int(sys.argv[1])
    #seed=1
    random.seed(seed)    
    np.random.seed(seed)
    
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
        
        #a = float(recovDistParams[0])
        #scale = float(recovDistParams[1])    
        
        return stats.gamma.pdf(u,a=a,scale=scale)
    
    
    def inf_distr(u,*CIdistParms):
         a = float(CIdistParms[0])
         scale = float(CIdistParms[1])
         tol = 1e-10
         #Basically: use de l'hopital when the ratio becomes 0/0
         #Otherwise go with definition. This regularises a lot the numerics
         
         x = np.where(stats.gamma.cdf(u,a=a,scale=scale)>1-tol,
                 1/scale - (a-1)/u,
                 stats.gamma.pdf(u,a=a,scale=scale)/(1- stats.gamma.cdf(u,a=a,scale=scale)))  
         x[0]=x[1]
         return x
    
    time_I_total = np.load('../Data/India_infected.npy')
    time_R_total = np.load('../Data/India_rec.npy')
    
    time_I_extended = np.sort(time_I_total)
    time_R_extended = np.sort(time_R_total)
    
    cases=len(time_I_total)
    grids=3000 
    result_x=[4.06914629e-04, 2.81346094e+00, 1.57537281e+00, 5.69340232e+00,
            1.89918160e+01]

    pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1], result_x[2]],\
                              recovDist=rec_haz, recovDistParms=[result_x[3],result_x[4]],\
                                  nTgrid=grids, nUgrid=grids, T=149)
    
    initialcondition1=np.exp(-pde.tgrids)
    
    X,Y, Yts=pde.finDiffUpdate(initialcond=initialcondition1, Yts=True)    
     
    ft = X*Yts/(1-X[-1])
    
    #Renormalize ft because its sum is slightly lower than 1
    ft = ft/(np.sum(ft)*pde.dx)
    
    
    infection_times=np.sort(np.random.choice(pde.tgrids, p=pde.dx*ft, size=3000))
    
    
    recovery_times=np.random.gamma(shape=result_x[3]**2/result_x[4], scale=result_x[4]/result_x[3], size=len(infection_times))

    recovery_times = infection_times + recovery_times 

    ll=log_likelihood_models(grids,hazard_inf=inf_distr,hazard_rec=rec_haz, rec_distr = rec_distr, 
                              T=135, infect_times=infection_times, recov_times=recovery_times , hazard_inf_par=2,rec_parms=2)
    
    #ll=log_likelihood_models(inf_distr,grids,hazard_rec=rec_haz, rec_distr = rec_distr, 
    #                          T=day[-1], infect_times=time_I_extended, recov_times=time_R_extended)
    
    result = ll.minimize_likelihood(np.array([1e-4,1,1,3,2]), np.array([1e-2,10,6,12, 25]), maxiter=20,swarmsize=100)

    np.savetxt("../Bootstrap/likelihood_%d.txt"%seed,result.x,newline=" ")
