#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 09:15:37 2021

@author: fra
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import sys
sys.path.append('../')

from Likelihood import log_likelihood_models
from PDE_solver import SIR_PDEroutine
from Selkealgo import Sellke_algo


#This is a fully workout example that uses all the methods of this repo

#It is divided into three parts:
    
# 1) Generate some Data according to some parameters
# 2) Solve the PDE with the same parameters to compare it to data
# 3) Use the likelihood to infer the parameters from the data and plot them

#These three bits are independent, so if one is interested only in solving the 
#PDE, they can copy and adapt the relative bit on a new code

if __name__ == "__main__":        

        np.random.seed(3) #Initialize random seed for reproducibility

        #1) Generate some data. In this example, we choose gamma distribution for
        #infectious period and exponential distribution for contact intervals
        N=10000
        T_f = 150   
        
        beta = 0.2          #infectiousness
        
        mean = 9             #mean infectious period
        variance =6          #variance of infectious period distribution
        
        scale = variance/mean  #inverserate
        a = mean**2/variance   #shape
        

        I_0 = 20              #Initial number of infected people
        tau = beta/(N-1)      #This is because there is a factor N in the Sellke construction
        check_data=False
        
        #First thing one needs is to generate the infectious periods for each node
        #Not all the nodes will use that, as not everyone will get infected.
        T_recover = stats.gamma.rvs(a=a,scale=scale, size=N)

        while check_data==False:
            time,I,S,R=Sellke_algo(tau,I_0,N,T_f,T_recover,showplot=False,return_full=False)
            if len(time)>200: #Make sure this epidemic is not dying out
                check_data=True

        plt.figure()
        plt.plot(time,I/N, color='black', label='data') #Plot Data
        plt.xlim(0,150)
        plt.ylim(0)
        #If you want to save the data in a Dataframe format
        #time, I, S
        
        #data= np.c_[np.array(time),np.array(I),np.array(S)]
        #name = "Example data,csv"
        #np.savetxt(name,zipped, header='time,I,S')        
        
        #======================================================================
        #2) solution of the PDE with the true parameters

        #We need two quantities: 
            #1) infectious period/recovery time hazard function
            #2) contact interval hazard functions
            
            #Note, in general one can use the fact that the hazard function is
            # pdf/survival function
            
        def rec_haz(u, *recovDistParams):    
            a = float(recovDistParams[0])
            scale = float(recovDistParams[1])  
            tol = 1e-10
            #Basically: use de l'hopital when the ratio becomes 0/0
            #Otherwise go with definition. This regularises a lot the numerics
            
            x = np.where(stats.gamma.cdf(u,a=a,scale=scale)>1-tol,
                         1/scale - (a-1)/u,
                         stats.gamma.pdf(u,a=a,scale=scale)/(1- stats.gamma.cdf(u,a=a,scale=scale)))
            return x
      
        def inf_haz(u,*CIdistParms):
            beta = float(CIdistParms[0])
            return beta*np.ones_like(u)     
        
        grids=T_f*20 #How fine the time solver grid (finder -> more precise but more time)
        
        rho = I_0/(N-I_0)
        
        
        pde = SIR_PDEroutine(rho, CIdist=inf_haz, CIdistParms=[beta],\
                                  recovDist=rec_haz, recovDistParms=[a, scale],\
                                      nTgrid=grids, nUgrid=grids, T=T_f)
        #Initial condition is a vector as long as the grid that contains the distribution
        #of recovery times of initially infected individuals 
        #In this case should be a delta in 0.
        initial_condition=np.zeros_like(pde.tgrids)
        initial_condition[0]=1
        #Solve the PDE
        S_pde,I_pde=pde.finDiffUpdate(intiial_condition=initial_condition)        
        
        plt.plot(pde.tgrids,I_pde, color='b', label= 'PDE')
    
        
        #======================================================================
        
        
        #3) Maximise the likelihood with infection and recovery times
        #We use infection and recovery times from the data generated
        infection_times=time[np.where(np.diff(I)>0)]
        recovery_times=time[np.where(np.diff(I)<0)]
        
        #We need also the recovery distribution to run the likelihood
        
        def rec_distr(u, *recovDistParams):
           a = float(recovDistParams[0])
           scale = float(recovDistParams[1])    
           
           return stats.gamma.pdf(u,a=a,scale=scale)        
        
        
        ll=log_likelihood_models(grids,hazard_inf=inf_haz,hazard_rec=rec_haz, 
                                 rec_distr = rec_distr, 
                              T=T_f, infect_times=infection_times,recov_times=recovery_times,hazard_inf_par=1,rec_parms=2)
        result = ll.minimize_likelihood(np.array([5e-4,0.01,1,0.1]), np.array([1e-2,2,20,1]))
                
        parameters=result.x


        #Plot the MLE
        
        pde = SIR_PDEroutine(parameters[0], CIdist=inf_haz, CIdistParms=[parameters[1]],\
                                  recovDist=rec_haz, recovDistParms=[parameters[2], parameters[3]],\
                                      nTgrid=grids, nUgrid=grids, T=T_f)
        #Initial condition in this case should be a delta in 0.
        initial_condition=np.zeros_like(pde.tgrids)
        initial_condition[0]=1
        #Solve the PDE
        S_mle,I_mle=pde.finDiffUpdate(intiial_condition=initial_condition)                    

        plt.plot(pde.tgrids,I_mle, color='r', label= 'MLE')
        
        plt.legend()
        