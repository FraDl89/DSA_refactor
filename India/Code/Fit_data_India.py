#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:24:28 2021

@author: fra
"""


import numpy as np
import sys
import scipy.stats as stats
import random

import pandas as pd
import matplotlib.pyplot as plt

sys.path.append('../../')

from Likelihood import log_likelihood_models
from PDE_solver import SIR_PDEroutine



if __name__=="__main__":
    np.random.seed(20131989)
    random.seed(4011994)
    
    
    
    path='../Data/case_time_series.csv'
    full_data = pd.read_csv(path)
    date_start = '2021-02-15' #first day to consider
    date_end = '2021-07-01' #last day (excluded)
    
    #date_start = '2020-03-01'
    #date_end = '2020-10-20'
    
    
    day_start= np.where(full_data["Date_YMD"] == date_start)[0][0]
    day_end = np.where(full_data["Date_YMD"] == date_end)[0][0]
    #These instructions are specific to the data set, should be self-explanatory 
    date_datetime = np.array(full_data["Date_YMD"][day_start:day_end])
    date_days = np.arange(day_start,day_end,1)
    
    infected = np.array(full_data["Daily Confirmed"][day_start:day_end])
    recov_deceased = np.array(full_data["Daily Recovered"][day_start:day_end])+ \
                    np.array(full_data["Daily Deceased"][day_start:day_end])
    total_recov = np.array(full_data["Total Recovered"][day_start:day_end])+ \
                    np.array(full_data["Total Deceased"][day_start:day_end])
    total_deaths = np.array(full_data["Total Deceased"][day_start:day_end])
    total_infected = np.array(full_data["Total Confirmed"][day_start:day_end])
    deaths = np.array(full_data["Daily Deceased"][day_start:day_end])
    
    recov =  np.array(full_data["Daily Recovered"][day_start:day_end])
    tot_rec =  np.array(full_data["Total Recovered"][day_start:day_end])
    
    
    dataframe = pd.DataFrame({'I':infected, 'I_t':total_infected, \
                              'RD':recov_deceased,'RD_t':total_recov,\
                            'R':recov,'R_t':tot_rec, 'D':deaths, \
                                'D_t':total_deaths},index=date_datetime)
        
    
    day = np.arange(day_start, day_end,1)
    day = day - day[0]  #Shift of initial day to match PDE
    
    #Distribute infection and recovery times
    #uniformly through the days (to avoid having atoms of probability)
    time_I_extended = np.zeros(sum(infected))
    time_R_extended = np.zeros(sum(recov_deceased))
    position=0
    '''
    
    
    #Room for improvement here, just a quick and dirty loop

    for dailycases in range(len(infected)):
        random_times = np.sort(np.random.uniform(0,1, infected[dailycases]))
        for time in random_times:
            if dailycases ==0:
                time_I_extended[position] = time
            else:
                time_I_extended[position] = day[dailycases-1]+time*(day[dailycases]-day[dailycases-1])
            position = position+1
    
    position=0
    
    for dailycases in range(len(recov_deceased)):
        random_times = np.sort(np.random.uniform(0,1, recov_deceased[dailycases]))
        for time in random_times:
            if dailycases ==0:
                time_R_extended[position] = time
            else:
                time_R_extended[position] = day[dailycases-1]+time*(day[dailycases]-day[dailycases-1])
            position = position+1
    
    np.save("Data/India_infected",time_I_extended)
    np.save("Data/India_rec",time_R_extended)
    '''
    time_I_total = np.load('../Data/India_infected.npy')
    time_R_total = np.load('../Data/India_rec.npy')
    
    time_I_extended = np.sort(time_I_total)
    time_R_extended = np.sort(time_R_total)
    
    
    
    
    #take 10000 samples equispaced in the space of indices
    size = 3000
    time_I_extended = np.sort(np.random.choice(time_I_extended,size,replace=False))
    
    time_R_extended = np.sort(np.random.choice(time_R_extended,size,replace=False))
    
    #sample_every = int(time_I_total.size/size)
    #time_I_extended = time_I_extended[::sample_every]
    #sample_every = int(time_R_total.size/size)
    #time_R_extended = time_R_extended[::sample_every]
    
    
    
    
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
    
    
    
    grids=3000
    
    ll=log_likelihood_models(grids,hazard_inf=inf_distr,hazard_rec=rec_haz, rec_distr = rec_distr, 
                              T=day[-1], infect_times=time_I_extended, recov_times=time_R_extended, hazard_inf_par=2,rec_parms=2)
    

    result = ll.minimize_likelihood(np.array([1e-4,1,1,3,2]), np.array([1e-2,10,6,12, 25]), maxiter=100,swarmsize=500)
    
    #print(result)
    
    print(result)
    #result_x=[4.06914629e-04, 2.81346094e+00, 1.57537281e+00, 5.69340232e+00,
    #        1.89918160e+01]
    # 26642.522743303794)
    result_x=result[0]
    
    pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1], result_x[2]],\
                              recovDist=rec_haz, recovDistParms=[result_x[3],result_x[4]],\
                                  nTgrid=grids, nUgrid=grids, T=day[-1])
    
    initialcondition1=np.exp(-pde.tgrids)
    
    X,Y=pde.finDiffUpdate(initialcond=initialcondition1)
        
    
    '''
    #Plots: I(t)
    plt.figure()
    plt.plot(pde.tgrids,Y, color='blue')
    plt.xlabel('time')
    plt.ylabel('cases')   
    '''
    #Plots: Recovery distribution    
    plt.figure()
    x = np.linspace(0,30,grids)
    plt.plot(x,rec_distr(x,*[result_x[3],result_x[4]]), label = 'Gamma: mean %.3f, var %.3f)'%(result_x[3],result_x[4]))
    plt.title("Foot and mouth")
    plt.xlabel('time')
    plt.legend()
    plt.ylabel('Recovery pdf')    
    plt.title('time to recovery')
    
    #Plots: Infectiousness distribution
    x = np.linspace(0,30,grids)
    
    plt.figure()
    plt.plot(x,stats.gamma.pdf(x,a=result_x[1],scale=result_x[2]), label = 'Gamma: mean %.3f, var %.3f'%(result_x[1]*result_x[2],result_x[1]*result_x[2]**2))
    plt.xlabel('time')
    plt.ylabel('Infection pdf')
    plt.legend()
    plt.title("infectiousness in time")
    #Plots: f_t
    from scipy.stats import gaussian_kde
    dense = gaussian_kde(time_I_extended)
    denseval = list(dense(x)  for x in day)        
    x=np.linspace(0,day,3000)
    plt.figure()
    
    empirical_density = -(np.diff(X)/pde.dx) /(1-X[-1])
    plt.plot(pde.tgrids[1:], empirical_density, color='b',label='MLE density uniform in cond' )
    plt.scatter(day,denseval, color='black', label = 'Empirical density')
    plt.xlabel("T")
    plt.ylabel("Conditional density")
    plt.legend()   
    
    
    
    #Gamma gamma
    #Likelihood 26557
    #[4.93174744e-04, 5.05115136e+00, 1.00000000e+00, 6.05378065e+00,1.93962008e+01]
    
    
    #Likelihood 26505
    
    #result_x=[5.44645257e-04, 3.59625285e+00, 1.64603099e+00, 8.21525320e+00,2.11575966e+01]
    
    #These are the best results for the first 250 days
    #GAMMA
    #array([3.16410516e-05, 1.41516534e-01, 8.36858743e+00, 7.32729328e+00])
    #date_start = '2020-03-01'
    #date_end = '2020-10-20'
    
    #[2.75066276e-04, 2.07164788e-01, 6.42516240e+00, 1.05904792e+01]
    #date_start = '2021-02-15' #first day to consider
    #date_end = '2021-06-15' #last day (excluded)
        
    #array([8.63431947e-04, 2.10000000e+00, 6.83694671e+00, 7.82871068e+00,
    #       1.58827459e+01])
    
    #[5.92080750e-04 1.50000000e+00 5.60343470e+00 6.91161839e+00
    # 2.10000000e+01] 
    
    
    '''
    import matplotlib.pyplot as plt
    pdeObj_1 = SIR_PDEroutine(result[0], CIdist=inf_distr, CIdistParms=[result[1]],\
                              recovDist=rec_haz, recovDistParms=[result[2],result[3]],\
                                  nTgrid=10000, nUgrid=10000, T=day[-1])
    X,Y=pdeObj_1.finDiffUpdate()
    
    #Effective N for DSA
    effective_N=(total_deaths[-1]+total_infected[-1])/(1-X[-1])
    print(effective_N/10**6)
    
    plt.plot(day, total_deaths+total_infected, color='r', label='cumulative cases observed')
    plt.plot(pdeObj_1.tgrids, effective_N*(1-X), label='estimate from PDE')
    plt.xlabel("Time (days)")
    plt.ylabel("Cumulative cases")
    plt.legend()
    #plt.savefig("India_fit_pde.pdf", format='pdf')
    plt.figure()
    #This is - \dot{S}/(1-X[-1]), should give the empirical density of infected
    empirical_density = -(np.diff(X)/pdeObj_1.dx) /(1-X[-1])
    
    from scipy.stats import gaussian_kde
    dense = gaussian_kde(time_I_extended)
    denseval = list(dense(x)  for x in day)
    
    plt.plot(pdeObj_1.tgrids[1:], empirical_density, color='b',label='PDE density' )
    plt.scatter(day,denseval, color='r', label = 'Empirical density')
    plt.xlabel("T")
    plt.ylabel("Conditional density")
    plt.legend()
    
    #Plot the recovery distribution as inferred from ML
    plt.figure()
    x = np.linspace(0,20,1000)
    plt.plot(x,rec_distr(x,*[result[2],result[3]]), label = 'Weibull(%.3f, %.3f)'%(result[2],result[3]))
    plt.title("India - first 250 days")
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('pdf')
    plt.xlim(0,30)
    '''

