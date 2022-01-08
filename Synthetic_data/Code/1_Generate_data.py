#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 09:15:37 2021

@author: fra
"""
import numpy as np
import gc
import pickle
import scipy.stats as stats
import scipy.special
import sys
from Selkealgo import Sellke_algo


#This bit generates synthetic datasets that are needed to generate estimates
#for section  4.1



if __name__ == "__main__":        

        seed = int(sys.argv[1])
        outputdir = sys.argv[2]
        #seed = 817
        #outputdir=''
        T_f = 80   

    
        np.random.seed(seed) #Initialize random seed for reproducibility

            
        N=10000
        I_0 = 50
        
        beta = 0.25          #infectiousness
        
        mean = 9             #mean infectious period
        variance =6          #variance of infectious period distribution
        
        scale = variance/mean  #inverserate
        a = mean**2/variance   #shape
        
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
            
        def rec_distr(u, *recovDistParams):
           a = float(recovDistParams[0])
           scale = float(recovDistParams[1])    
           
           return stats.gamma.pdf(u,a=a,scale=scale)
       
        def inf_distr(u,*CIdistParms):
            beta = float(CIdistParms[0])
            return beta*np.ones_like(u)     
            

        I_0 = 50
        tau = beta/(N-1)
        check_data=False
        T_recover = stats.gamma.rvs(a=a,scale=scale, size=N)

        while check_data==False:
            time,I,S,R,times_to_sample,grid_I,grid_S=\
            Sellke_algo(tau,I_0,N,T_f,T_recover,showplot=False,return_full=True)
            if len(time)>200: #Make sure this epidemic is not dying out
                check_data=True


        #import matplotlib.pyplot as plt
        #plt.plot(time,I)
        
        zipped = (np.array(time),np.array(I),np.array(S))
        name = "%s/Gamma_distr_seed%d"%(outputdir,seed)
        np.save(name,zipped)        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

                
