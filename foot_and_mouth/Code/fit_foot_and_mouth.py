#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2020-2022 Francesco Di Lauro. All Rights Reserved.
See Licence file for details.
"""

import numpy as np
import matplotlib.pyplot as plt
#from untitled0 import log_likelihood_models
import scipy.stats as stats
import random
from pyDOE import lhs

import sys
sys.path.append('../../')

from Likelihood import log_likelihood_models
from PDE_solver import SIR_PDEroutine


'''
This code is used to fit the real dataset from the foot and mouth epidemics.
The processed_time_daily_cases.txt might be available upon reasonable request
'''

if __name__=="__main__":
     
    np.random.seed(1404) #Initialize random seed for reproducibility
    processed_time_cases = '../Data/processed_time_daily_cases.txt'  #Folder with FMD data
    
    #Cases we are interested in
    #plt.figure()        
    day,cases = np.loadtxt(processed_time_cases, delimiter='  ',unpack=True)    
    
    #plt.scatter(day,cases, color='black', marker='x')

    #plt.xlabel(r'$t$ (day)')
    #plt.ylabel(r'cases')
    
    
    day = day.astype(int)
    
    T_f = 40
    day=day[:T_f]
    cases = cases[:T_f]


    #plt.scatter(day,cases, color='red', marker='x')
    
    day = day-day[0] #make time begin from day 0
    
    time_extended = np.zeros(sum(cases.astype(int)))
    
    position=1
    #Distribute uniformly through the days the times of infections
    for dailycases in range(1,len(cases)):
        random_times = np.sort(np.random.uniform(0,1, np.int(cases[dailycases])))
        for time in random_times:
            time_extended[position] = day[dailycases-1]+time*(day[dailycases]-day[dailycases-1])
            position +=1
    
    #N = 131000 #Approximately the number of farms in the UK as per Ferguson's paper
               #https://pubmed.ncbi.nlm.nih.gov/11303090/
               #The Foot-and-Mouth Epidemic in Great Britain: Pattern of Spread 
               #and Impact of Interventions http://science.sciencemag.org/content/292/5519/1155
    
 
    #Gamma
    '''
    def rec_haz(u, *recovDistParams):    
         a = float(recovDistParams[0])**2/float(recovDistParams[1])
         scale = float(recovDistParams[1])/float(recovDistParams[0]) 
         tol = 1e-10
         #Basically: use de l'hopital when the ratio becomes 0/0
         #Otherwise go with definition. This regularises a lot the numerics
         x=np.zeros_like(u)
         bad_indexes=stats.gamma.cdf(u,a=a,scale=scale)>1-tol
         x[bad_indexes]=1/scale - (a-1)/u[bad_indexes]
         x[~bad_indexes]= stats.gamma.pdf(u[~bad_indexes],a=a,scale=scale)/(1- stats.gamma.cdf(u[~bad_indexes],a=a,scale=scale))

 
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
   
    '''
    def rec_haz(u,*recovDistParams):
        rate=float(recovDistParams[0])
        return rate*np.ones_like(u)

    def rec_distr(u,*recovDistParams):
        rate=float(recovDistParams[0])
        return stats.expon.pdf(u,scale=1/rate)
    def inf_distr(u,*CIdistParms):
        rate=float(CIdistParms[0])
        return rate*np.ones_like(u)
        

    
    grids=1600
    ll=log_likelihood_models(grids,hazard_inf=inf_distr,hazard_rec=rec_haz, rec_distr = rec_distr, 
                              T=T_f, infect_times=time_extended, hazard_inf_par=1,rec_parms=1)

    #result = ll.minimize_likelihood(np.array([5e-4,2.1,2,2,1]), np.array([1e-2,10,10,9, 21]), maxiter=200,swarmsize=150)
    result = ll.minimize_likelihood(np.array([1e-4,0.1,0.05]),np.array([1e-2,0.9,0.9]),maxiter=200,swarmsize=150)
    print(result)
    #result expected:
    #result_x=[8.23024809e-03, 2.13623063e+00, 4.75098558e+00, 4.97839683e+00,1.08327439e+01]
    plt.figure()
    
    result_x=result[0]


        

    #pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1], result_x[2]],\
    #                          recovDist=rec_haz, recovDistParms=[result_x[3],result_x[4]],\
    #                              nTgrid=grids, nUgrid=grids, T=T_f)
    pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1]],\
                             recovDist=rec_haz, recovDistParms=[result_x[2]],\
                                 nTgrid=5000, nUgrid=5000, T=T_f)    
    initialcondition1=np.exp(-pde.tgrids)
    
    X,Y=pde.finDiffUpdate(initialcond=initialcondition1)
 
    N_eff = sum(cases)/(1-X[-1])
    rho_0 = result_x[0]*N_eff
    
    #Plots: I(t)
    '''
    plt.plot(pde.tgrids,Y, color='blue')
    plt.xlabel('time')
    plt.ylabel('cases')   
    
    #Plots: Recovery distribution    
    plt.figure()
    x = np.linspace(0,24,grids)
    plt.plot(x,rec_distr(x,*[result_x[3],result_x[4]]), label = 'Gamma: mean %.3f, var %.3f)'%(result_x[3],result_x[4]))
    #plt.plot(x,rec_distr(x,*[result_x[2]]), label = 'exponential: rate %.2f'%(result_x[2]))
    plt.title("Foot and mouth")
    plt.xlabel('time')
    plt.legend()
    plt.ylabel('Recovery pdf')    

    #Plots: Infectiousness distribution
    
    plt.figure()
    plt.plot(x,stats.weibull_min.pdf(x,result_x[1], scale=result_x[2]), label = 'Weib(%.3f, %.3f)'%(result_x[1],result_x[2]))
    #plt.plot(x,stats.expon.pdf(x,scale=1/result_x[1]),label = 'exponential: rate %.2f'%(result_x[1]))
    plt.xlabel('time')
    plt.ylabel('Infection pdf')
    '''
 
    #Plots: f_t
    from scipy.stats import gaussian_kde
    dense = gaussian_kde(time_extended)
    denseval = list(dense(x)  for x in day)        
    x=np.linspace(0,day,3000)
    plt.figure()

    empirical_density = -(np.diff(X)/pde.dx) /(1-X[-1])
    plt.plot(pde.tgrids[1:], empirical_density, color='b',label='MLE density uniform in cond' )
    plt.scatter(day,denseval, color='black', label = 'Empirical density')
    plt.xlabel("T")
    plt.ylabel("Conditional density")
    plt.legend()    
    
    #N_eff and initial number of infected
    N_eff = sum(cases)/(1-X[-1])
    rho_0 = result_x[0]*N_eff

    np.savetxt('Results/estimates_FMD_exp_exp_T_f_%d.txt'%T_f, np.c_[result_x[0],result_x[1],result_x[2],rho_0,N_eff, sum(cases)])
    #np.savetxt('Results/estimates_FMD_Weib_Gamma_T_f_%d.txt'%T_f, np.c_[result_x[0],result_x[1],result_x[2],result_x[3],result_x[4],rho_0,N_eff, sum(cases)])    
    
    #np.savetxt('fit_empirical_density_foot_and_mouth.txt', np.c_(pde.tgrids[1:],empirical_density))
    #np.savetxt('empirical_density_foot_and_mouth.txt', np.c_(day,denseval))
    
  
    
    '''

    Unidentifiability:
        
    result_x_gam_weib=[7.49589223e-03, 3.25712686e+00, 4.70202301e+00, 2.96796282e+00,
           9.87841682e+00]
    
    result_x_2_gam_weib=[5.19105809e-03, 1.00027262e+00, 3.94524974e+00, 5.38972307e+00,
           8.99641505e+00]
    
    #These two sets of data have roughly the same likelihood but look very different



    result_x = result_x_gam_weib

    pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1], result_x[2],result_x[3]],\
                              recovDist=rec_haz, recovDistParms=[result_x[4],result_x[5]],\
                                  nTgrid=grids, nUgrid=grids, T=day[-1])
    
    initialcondition0='Unif'
   
    initialcondition1=np.exp(-pde.tgrids)
    
    X,Y=pde.finDiffUpdate(initialcond=initialcondition)
    plt.plot(pde.tgrids,Y, color='blue')
    plt.xlabel('time')
    plt.ylabel('cases')
    
    X,Y=pde.finDiffUpdate(initialcond=initialcondition1)
    plt.plot(pde.tgrids,Y, color='red')
        
    
    plt.figure()
    x = np.linspace(0,24,grids)
    plt.plot(x,rec_distr(x,*[result_x[4],result_x[5]]), label = 'Gamma: mean %.3f, var %.3f)'%(result_x[4],result_x[5]))
    plt.title("Foot and mouth")
    plt.xlabel('time')
    plt.legend()
    plt.ylabel('Recovery pdf')
    plt.figure()
    plt.plot(x,stats.gamma.pdf(x,result_x[1], scale=result_x[2], loc=result_x[3]), label = 'Weib(%.3f, %.3f)'%(result_x[1],result_x[2]))
    plt.xlabel('time')
    plt.ylabel('Infection pdf')
  
    
    from scipy.stats import gaussian_kde
    dense = gaussian_kde(time_extended)
    denseval = list(dense(x)  for x in day)        
    x=np.linspace(0,day,3000)
    plt.figure()

    empirical_density = -(np.diff(X)/pde.dx) /(1-X[-1])
    plt.plot(pde.tgrids[1:], empirical_density, color='b',label='MLE density uniform in cond' )
    plt.scatter(day,denseval, color='r', label = 'Empirical density')
    plt.xlabel("T")
    plt.ylabel("Conditional density")
    plt.legend()    
    
    
    #survival function
    ft = (X-X[-1])/(1-X[-1])
    plt.plot(pde.tgrids,ft)
    
    N_eff = sum(cases)/(1-X[-1])
    rho_0 = result_x[0]*N_eff
    
    #Solve in the future
    pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1]],\
                              recovDist=rec_haz, recovDistParms=[result_x[2],result_x[3]],\
                                  nTgrid=grids, nUgrid=grids, T=120)
    X,Y=pde.finDiffUpdate()    
    plt.plot(pde.tgrids,X*(1-rho_0/N_eff)*N_eff)
    '''
    
    
