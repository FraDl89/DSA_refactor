#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:30:09 2021

@author: fra
"""



import numpy as np
from Likelihood import log_likelihood_models
import scipy.stats as stats
from PDE_solver import SIR_PDEroutine
import random
import sys
if __name__=="__main__":
     
    seed =  int(sys.argv[1])
    
    random.seed(seed)
    np.random.seed(seed)
    processed_time_cases = 'foot_and_mouth/processed_time_daily_cases.txt'  #Folder with FMD data
    
        
    day,cases = np.loadtxt(processed_time_cases, delimiter='  ',unpack=True)    
    
    day = day.astype(int)
    
    T_f = 80
    day=day[:T_f]
    cases = cases[:T_f]
    #number of cases
    tot_cases=int(sum(cases))
    day = day-day[0] #make time begin from day 0
    

    #Gamma
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
           #Weibull
           shape=CIdistParms[0]
           scale=CIdistParms[1]
           return shape/scale*(u/scale)**(shape-1)
       
    #This is the output from the model     
    result_x=[8.23024809e-03, 2.13623063e+00, 4.75098558e+00, 4.97839683e+00,1.08327439e+01]
    grids=1600
        
    total_extractions=int(sum(cases))
       

    pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1], result_x[2]],\
                              recovDist=rec_haz, recovDistParms=[result_x[3],result_x[4]],\
                                  nTgrid=grids, nUgrid=grids, T=day[-1])
    
    initialcondition1=np.exp(-pde.tgrids)
    
    X,Y, Yts=pde.finDiffUpdate(initialcond=initialcondition1, Yts=True)


    ft = X*Yts/(1-X[-1])
    
    #Renormalize ft because its sum is slightly lower than 1
    ft = ft/(np.sum(ft)*pde.dx)
    
    
    infection_times=np.sort(np.random.choice(pde.tgrids, p=pde.dx*ft, size=tot_cases))
    
    
    ll=log_likelihood_models(grids,hazard_inf=inf_distr,hazard_rec=rec_haz, rec_distr = rec_distr, 
                              T=T_f, infect_times=infection_times, hazard_inf_par=2,rec_parms=2)

    result = ll.minimize_likelihood(np.array([5e-4,2.1,2,2,1]), np.array([1e-2,10,10,9, 21]))
    print(*result.x)


        