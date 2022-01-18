#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2020-2022 Francesco Di Lauro. All Rights Reserved.
See LICENCE file for details
"""

import scipy.stats as stats
import random
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import convolve

sys.path.append('../../')

from PDE_solver import SIR_PDEroutine
from Likelihood import log_likelihood_models

if __name__=="__main__":

    
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




    path='../Data/case_time_series.csv'
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


    

    
    
    
    big_matrix_I=np.zeros((500,3000-1))
    
    big_matrix_R=np.zeros((500,3000-1))
 
    #matrix_ij refers to the bootstrap extraction i, time j of the pde
       
    R_0distr=np.zeros(500)
    
    for n in range(1,501):

            
        result_x=np.loadtxt("../Bootstrap/likelihood_%d.txt"%n)

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


        x=np.linspace(0, 50, 2000)
        a = result_x[3]**2/result_x[4]
        scale = result_x[4]/result_x[3]        
        
        prob_infectious=1-stats.gamma.cdf(x,a=a,scale=scale)
        rate_transmission=stats.gamma.pdf(x,a=result_x[1], scale=result_x[2])/(1- stats.gamma.cdf(x,a=result_x[1], scale=result_x[2]))       
        R_0distr[n-1]= np.sum(rate_transmission*prob_infectious*50/2000)
       
        
       
    two_five_I,fifty_I,nine_five_I = np.quantile(big_matrix_I, [0.025,0.5,0.975], axis=0)

    two_five_R,fifty_R,nine_five_R = np.quantile(big_matrix_R, [0.025,0.5,0.975], axis=0)

    np.savetxt("../Data/infected_confidence_intervals.txt",np.c_[pde.tgrids[1:],two_five_I,fifty_I,nine_five_I])
    np.savetxt("../Data/recovered_confidence_intervals.txt",np.c_[pde.tgrids[1:],two_five_R,fifty_R,nine_five_R])
    
    np.savetxt("../Data/infected_cases.txt",np.c_[day,infected])
    np.savetxt("../Data/recovered_recoveries.txt",np.c_[day,recov_deceased])

    np.savetxt("../Data/R_0_distr.txt", np.c_[R_0distr])
    
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

    '''     
    
