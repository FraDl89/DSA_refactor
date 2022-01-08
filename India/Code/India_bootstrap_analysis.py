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
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import convolve

if __name__=="__main__":
    def moving_average(a,n):
        N=len(a)
        return np.array([np.mean(a[i:i+n]) for i in np.arange(0,N-n+1)])
    
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
    
    time_I_total = np.load('India/India_infected.npy')
    time_R_total = np.load('India/India_rec.npy')
    
    time_I_extended = np.sort(time_I_total)
    time_R_extended = np.sort(time_R_total)
    
    cases=len(time_I_total)    




    path='India/case_time_series.csv'
    full_data = pd.read_csv(path)
    date_start = '2021-02-15' #first day to consider
    date_end = '2021-07-01' #last day (excluded)
    

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


    mov_avg_I=moving_average(infected,7)
    mov_avg_R=moving_average(recov_deceased, 7)
    
    differences_I=mov_avg_I-infected[6:]
    differences_R=mov_avg_R-recov_deceased[6:]    
    
    var_diff_I= np.var(differences_I)
    var_diff_R= np.var(differences_R)    
    
    
    
    big_matrix_I=np.zeros((500,3000-1))
    
    big_matrix_R=np.zeros((500,3000-1))
        
    for n in range(1,501):

            
        result_x=np.loadtxt("India/Bootstrap/likelihood_%d.txt"%n)

        grids=3000
        T=135
        pde= SIR_PDEroutine(result_x[0], CIdist=inf_distr, CIdistParms=[result_x[1], result_x[2]],\
                                  recovDist=rec_haz, recovDistParms=[result_x[3],result_x[4]],\
                                      nTgrid=grids, nUgrid=grids, T=T)
        
        rec_times=rec_distr(pde.tgrids,*[result_x[1],result_x[2]])     
        
        initialcondition1=np.exp(-pde.tgrids)
        
        
        X,Y, Yts=pde.finDiffUpdate(initialcond=initialcondition1, Yts=True)        
        
        #N_eff and initial number of infected
        N_eff = sum(infected)/(1-X[-1])
        rho_0 = result_x[0]*N_eff
         
        big_matrix_I[n-1]=-np.diff(X)/pde.dx *N_eff
        #plt.plot(pde.tgrids[1:],-np.diff(X)/pde.dx *N_eff, alpha=0.6)
        
        
        big_matrix_R[n-1] =  -(np.diff(Y)+ np.diff(X))/pde.dx *N_eff

       
        
       
    two_five_I,fifty_I,nine_five_I = np.quantile(big_matrix_I, [0.025,0.5,0.975], axis=0)

    two_five_R,fifty_R,nine_five_R = np.quantile(big_matrix_R, [0.025,0.5,0.975], axis=0)

    np.savetxt("India/infected_confidence_intervals.txt",np.c_[pde.tgrids[1:],two_five_I,fifty_I,nine_five_I])
    np.savetxt("India/recovered_confidence_intervals.txt",np.c_[pde.tgrids[1:],two_five_R,fifty_R,nine_five_R])
    
    np.savetxt("India/infected_cases.txt",np.c_[day,infected])
    np.savetxt("India/recovered_recoveries.txt",np.c_[day,recov_deceased])
    
    '''
    #Not moving average and not inflated conf int
    plt.figure()

    plt.scatter(day,infected,color='k', label='cases')
    
    plt.fill_between(pde.tgrids[1:], two_five_I, nine_five_I, color='b', alpha=0.4, label='95 conf interval') 
    plt.xlabel('day')
    plt.ylabel('cases')
    plt.xlim(0,80)
    plt.ylim(0,52)
    plt.legend()
    plt.title('foot and mouth')    
    
    
    plt.figure()

    plt.scatter(day,recov_deceased,color='k', label='recoveries')
    plt.fill_between(pde.tgrids[1:], two_five_R, nine_five_R, color='b', alpha=0.4, label='95 conf interval') 
    plt.xlabel('day')
    plt.ylabel('recoveries')
    plt.xlim(0,80)
    plt.ylim(0,52)

    
    plt.fill_between(pde.tgrids[1:], two_five_I, nine_five_I, color='b', alpha=0.4, label='95 conf interval') 
    plt.xlabel('day')
    plt.ylabel('cases')
    plt.xlim(0,80)
    plt.ylim(0,52)
    plt.legend()
    plt.title('foot and mouth')    
        
    
    

    #Moving average and not inflated conf int
    plt.figure()

    plt.scatter(day[6:],mov_avg_I,color='r', label='cases 7 day moving avg')
    
    plt.fill_between(pde.tgrids[1:], two_five_I, nine_five_I, color='b', alpha=0.4, label='95 conf interval') 
    plt.xlabel('day')
    plt.ylabel('cases')
    plt.xlim(0,80)
    plt.ylim(0,52)
    plt.legend()
    plt.title('foot and mouth')    
    

    
    
    
    #Moving average and inflated conf int
    plt.figure()
    centralpoint= (two_five_I + nine_five_I)/2
    distance= (nine_five_I-two_five_I)*np.sqrt(var_diff_I)
    plt.scatter(day[6:],mov_avg_I,color='k', label='7-day moving average')
    
    plt.fill_between(pde.tgrids[1:], centralpoint-0.5*distance, centralpoint+0.5*distance, color='b', alpha=0.4, label='95+ conf interval') 
    plt.xlabel('day')
    plt.ylabel('cases')
    plt.xlim(0,80)
    plt.ylim(0,52)
    plt.legend()
    plt.title('foot and mouth')
    
    plt.show()
    #np.savetxt("foot_and_mouth/inflted_confidence_intervals.txt",np.c_[pde.tgrids[1:],two_five,fifty,nine_five])
    '''     
    