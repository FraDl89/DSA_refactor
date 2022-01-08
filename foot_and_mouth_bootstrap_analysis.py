#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 14:06:27 2021

@author: fra
"""


import numpy as np
from Likelihood import log_likelihood_models
import scipy.stats as stats
from PDE_solver import SIR_PDEroutine
import random
import sys

import matplotlib.pyplot as plt

if __name__=="__main__":
    def moving_average(a,n):
        N=len(a)
        return np.array([np.mean(a[i:i+n]) for i in np.arange(0,N-n+1)])
    processed_time_cases = 'foot_and_mouth/processed_time_daily_cases.txt'  #Folder with FMD data
        
    day,cases = np.loadtxt(processed_time_cases, delimiter='  ',unpack=True)    
    
    
    T_f = 80
    day = day-day[0] #make time begin from day 0
    day_long=day
    day=day[6:T_f]

    cases = cases[:T_f]
    mov_avg=moving_average(cases,7)
    
    differences=mov_avg-cases[6:]
    var_diff= np.var(differences)
    
    
    
    
    big_matrix=np.zeros((500,1600-1))
    
    for n in range(1,501):
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
           
            
        with open('foot_and_mouth/%d.txt'%n) as f:
            values = f.readlines()[-1].strip().split()
            
            result_x = [float(i) for i in values]
        grids=1600
        T=80
        pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1], result_x[2]],\
                                  recovDist=rec_haz, recovDistParms=[result_x[3],result_x[4]],\
                                      nTgrid=grids, nUgrid=grids, T=T)
        
        
        
        initialcondition1=np.exp(-pde.tgrids)
        
        
        X,Y, Yts=pde.finDiffUpdate(initialcond=initialcondition1, Yts=True)        
        
        #N_eff and initial number of infected
        N_eff = sum(cases)/(1-X[-1])
        rho_0 = result_x[0]*N_eff
         
        big_matrix[n-1]=-np.diff(X)/pde.dx *N_eff
        #plt.plot(pde.tgrids[1:],-np.diff(X)/pde.dx *N_eff, alpha=0.6)
        
       
        
       
        
       
    two_five,fifty,nine_five = np.quantile(big_matrix, [0.025,0.5,0.975], axis=0)
        
    #Not moving average and not inflated conf int
    plt.figure()

    plt.scatter(day_long[:T_f],cases,color='k', label='cases')
    
    plt.fill_between(pde.tgrids[1:], two_five, nine_five, color='b', alpha=0.4, label='95 conf interval') 
    plt.xlabel('day')
    plt.ylabel('cases')
    plt.xlim(0,80)
    plt.ylim(0,52)
    plt.legend()
    plt.title('foot and mouth')    
    
    

    #Moving average and not inflated conf int
    plt.figure()

    plt.scatter(day,mov_avg,color='k', label='cases')
    
    plt.fill_between(pde.tgrids[1:], two_five, nine_five, color='b', alpha=0.4, label='95 conf interval') 
    plt.xlabel('day')
    plt.ylabel('cases')
    plt.xlim(0,80)
    plt.ylim(0,52)
    plt.legend()
    plt.title('foot and mouth')    
    

    
    
    
    #Moving average and inflated conf int
    plt.figure()
    centralpoint= (two_five + nine_five)/2
    distance= (nine_five-two_five)*np.sqrt(var_diff)
    plt.scatter(day,mov_avg,color='k', label='7-day moving average')
    
    plt.fill_between(pde.tgrids[1:], centralpoint-0.5*distance, centralpoint+0.5*distance, color='b', alpha=0.4, label='95+ conf interval') 
    plt.xlabel('day')
    plt.ylabel('cases')
    plt.xlim(0,80)
    plt.ylim(0,52)
    plt.legend()
    plt.title('foot and mouth')
    
    plt.show()
    #np.savetxt("foot_and_mouth/inflted_confidence_intervals.txt",np.c_[pde.tgrids[1:],two_five,fifty,nine_five])
    
    