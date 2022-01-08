#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:18:54 2020

@author: ld288
"""

import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
from heapq import *

'''
local implementation of the Sellke algorithm, following 
https://courses.helsinki.fi/sites/default/files/course-material/4672649/Sellke.pdf
'''


def Sellke_algo(tau,I_0,N,T_f,T_rec,showplot=True,return_full=True):
    '''
    args: 
        tau - scalar, per-contact infection rate
        I_0 - scalar, initial number of infected nodes
        N - scalar, total number of nodes
        T_f - scalar, final time 
        T_rec - np.array-like, distribution of recovery times
        showplot - boolean
        return_full - use True to return the timegrid as well - for averages
    '''
    #this is a grid used for the average
    times_to_sample = np.linspace(0,T_f,10*T_f)
    #used for the average
    average_sellke_S = np.zeros_like(times_to_sample)
    average_sellke_I = np.zeros_like(times_to_sample)
    

    #Generates N-I_0 recoveries
    T_infections = np.random.exponential(1,N-I_0)      #This is used for infections
    #sort infection times in ascending order
    T_infections = np.sort(T_infections)
    
    #We need to sort the first I_0 elements of T_infections
    #Generate the first I_0 recoveries    
    T_recover = T_rec.tolist() #This tells recoveries
    T_initial = T_rec[:I_0]
    T_recover = T_recover[I_0:]
    recqueue=[]
    for t in T_initial:
        #heappush to insert events in order
        heappush(recqueue,t)
    
    #total number of infected/susceptible
    I=[I_0]
    S = [N-I_0]
    #time of events
    time = [0]
    #this keeps track of the total pressure up to time t
    infpressuretotimet=0

    while(I[-1]>0 and time[-1]<T_f): #until no more events
        inf_pressure = tau*I[-1] #infectious pressure at time t
                
        #to find the next event at time t, we need to
        #find t such that totalpressure+ inf_pressure*delta_t = T_infections[0]
        #T_infections[0] is the first available infection time
        #or delta_t = (T_infections[0]-totalpressure)/inf_pressure
        try:
            delta_t_infection= (T_infections[0]-infpressuretotimet)/inf_pressure
        except:
            #if T_infections[0] does not exist -i.e. everyone is I or R
            delta_t_infection = np.inf
        #take the earliest recovery event and tell when will it happen from 
        #time t
        delta_t_first_recovery =  heappop(recqueue) - time[-1]
        #take the minimum of the two
        delta_t = min(delta_t_infection,delta_t_first_recovery)
        
        if  delta_t_infection < delta_t_first_recovery:
            #event is infection
            #now we know the time to the next event and we can update the total
            #infectous pressure accordingly
            infpressuretotimet +=delta_t_infection*I[-1]*tau
            #event is infection
            I.append(I[-1]+1) 
            S.append(S[-1]-1) 
            #put back in queue the recovery - we did not use it
            #note: unfortunately with heap there is no way of simply reading
            #the first event, you have to take it and then put it back. Weird
            heappush(recqueue,delta_t_first_recovery + time[-1]) 
            #update time
            time.append(time[-1]+delta_t)
            #This is very important, now we drop the first element of the 
            #infections vector, as we have used it
            T_infections=T_infections[1:]
            #take the first recovery event we generate and push it into
            #the ordered queue. Remember that recover[0] comes from a distribu-
            #tion that we take into account
            heappush(recqueue,T_recover[0]+time[-1]) 
            
            #shift recoveries by one
            T_recover = T_recover[1:]
        else: #event is recovery
            #similar as before
            infpressuretotimet +=delta_t*tau*I[-1]
            time.append(time[-1]+delta_t)
            #recovery
            I.append(I[-1]-1)
            S.append(S[-1])
    if showplot==True:
        #plot stuff
        plt.plot(time,I, color='gray')
        plt.plot(time,S, color='gray', linestyle='--') 
   

    #produce recovered by difference
    R = [N-I[t]-S[t] for t in range(len(time))]
    if return_full==True:
        for index,t_grid in enumerate(times_to_sample):
            #I take the closest event to time t, although I should be 
            #taking the earliest event closest to time t... it was simpler this
            #way, will modify eventually
            idx = (np.abs(time - t_grid )).argmin()
            average_sellke_I[index] += I[idx]
            average_sellke_S[index] += S[idx]  
        return time,I,S,R,times_to_sample,average_sellke_I,average_sellke_S
    else:
        return time,I,S,R










if __name__=="__main__":
    np.random.seed(2)

    N = 1000
    I_0 = 24
    #epidemic parameters
    tau = 0.15/N
    #Complete graph
    #G = nx.complete_graph(N)
    
    
    
    #this is a grid used for the average
    T_f = 40
    
    times_to_sample = np.linspace(0,T_f,10*T_f)
    
    '''
    Erlang plots
    '''
    #This is for erlang with low variance ~18
    average_erlang1_S = np.zeros_like(times_to_sample)
    average_erlang1_I = np.zeros_like(times_to_sample)
    #This is for erlang with high variance ~90    
    average_erlang2_S = np.zeros_like(times_to_sample)
    average_erlang2_I = np.zeros_like(times_to_sample)
    #this is for exponential with gamma 1/14 and variance ~192
    average_gamma_S = np.zeros_like(times_to_sample)
    average_gamma_I = np.zeros_like(times_to_sample)

    samples=50
    for i in range(samples):
        shape, scale = 10., 14/10.  # mean=14, std=10*14**2/100
        s_1 = np.random.gamma(shape, scale, 1000)
        shape, scale = 2, 14/2
        s_2 = np.random.gamma(shape, scale, 1000)    
        
        time,I,S,R,times_to_sample,grid_I,grid_S=\
        Sellke_algo(tau,I_0,N,T_f,s_1,showplot=False,return_full=True)
        average_erlang1_I += grid_I
        average_erlang1_S += grid_S    
#        if i==0:
#           plt.plot(time, S, color='gray',alpha=0.4, linestyle='--', label='Erlang low variance')    
#        else:
#            plt.plot(time,I, color='gray',alpha=0.4)
#            plt.plot(time,S, color='gray',alpha=0.4, linestyle='--')          
        time,I,S,R,times_to_sample,grid_I,grid_S=\
        Sellke_algo(tau,I_0,N,T_f,s_2,showplot=False,return_full=True)
        average_erlang2_I += grid_I
        average_erlang2_S += grid_S    

        gamma=1/14

        T_recover = np.random.exponential(1/gamma, N)#This tells recoveries
       
        time,I,S,R,times_to_sample,grid_I,grid_S=\
        Sellke_algo(tau,I_0,N,T_f,T_recover,showplot=False,return_full=True)
        average_gamma_I += grid_I
        average_gamma_S += grid_S    


#        if i==0:
#           plt.plot(time, S, color='r',alpha=0.4, linestyle=':', label='Erlang high variance')    
#        else:
#            plt.plot(time,I, color='r',alpha=0.4)
#            plt.plot(time,S, color='r',alpha=0.4, linestyle=':')          
        

    plt.plot(times_to_sample, average_erlang1_I/samples, color='k', label='average small-var erlang')
    plt.plot(times_to_sample, average_erlang1_S/samples, color='k')
    plt.plot(times_to_sample, average_erlang2_I/samples, color='b', label='average high-var erlang')
    plt.plot(times_to_sample, average_erlang2_S/samples, color='b')
    plt.plot(times_to_sample, average_gamma_I/samples, color='r', label='average exponential')
    plt.plot(times_to_sample, average_gamma_S/samples, color='r')
    plt.legend(loc='upper right')
    plt.xlabel(r'$t$')
    plt.ylabel(r'Prevalence')
    plt.xlim(0,T_f)
    plt.ylim(0,1000)
    print(np.mean(s_1),np.mean(s_2),np.mean(T_recover))
    print(np.var(s_1),np.var(s_2),np.var(T_recover))
    '''
     fast_SIR from Joel's EoN vs Sellke
    #used for the average
    average_sellke_S = np.zeros_like(times_to_sample)
    average_sellke_I = np.zeros_like(times_to_sample)
    tau = 2/N
    G = nx.complete_graph(N)
    average_joel_S = np.zeros_like(times_to_sample)
    average_joel_I = np.zeros_like(times_to_sample)
        
    
    gamma=1.8
    tau = 2.5/N
    #samples 
    samples=50
    for i in range(samples):
        T_recover = np.random.exponential(1/gamma, N)#This tells recoveries
       
        time,I,S,R,times_to_sample,grid_I,grid_S=\
        Sellke_algo(tau,I_0,N,T_f,T_recover,showplot=False,return_full=True)
        average_sellke_I += grid_I
        average_sellke_S += grid_S    
        #if i==0:
        #   plt.plot(time, S, color='gray', linestyle='--', label='Sellke SIR')    
        #else:
        #    plt.plot(time,I, color='gray')
        #    plt.plot(time,S, color='gray', linestyle='--') 

        #compare with EON
        initial_size = I_0
        
        t, S, I, R = EoN.fast_SIR(G, tau, gamma,
                                    initial_infecteds = range(initial_size))
        #if i==0:
        #   plt.plot(t, S, color='r', linestyle='--', label='Joel SIR')
        #else:
        #    plt.plot(t, I, color='r')
        #    plt.plot(t, S, color='r', linestyle='--')
        for index,t_grid in enumerate(times_to_sample):
            idx = (np.abs(t - t_grid )).argmin()
            average_joel_I[index] += I[idx]
            average_joel_S[index] += S[idx]
    
    plt.plot(times_to_sample, average_sellke_I/(samples*N), color='k', label='average sellke')
    plt.plot(times_to_sample, average_sellke_S/(samples*N), color='k')
    plt.plot(times_to_sample, average_joel_I/(samples*N), color='b', label='average joel')
    plt.plot(times_to_sample, average_joel_S/(samples*N), color='b')
    plt.legend(loc='upper right')
    plt.xlabel(r'$t$')
    plt.ylabel(r'Prevalence')
    
    
    
    np.random.seed(2)

    N = 1000
    I_0 = 2
    #epidemic parameters
    tau = 0.15/N
    #Complete graph
    
    #this is a grid used for the average
    T_f = 200
    
    times_to_sample = np.linspace(0,T_f,10*T_f)
    '''
    '''
    Erlang plots
    
    #This is for erlang with low variance ~18
    average_erlang1_S = np.zeros_like(times_to_sample)
    average_erlang1_I = np.zeros_like(times_to_sample)

    #samples=50
    #for i in range(samples):
    shape, scale = 10., 14/10.  # mean=14, std=10*14**2/100
    s_1 = np.random.gamma(shape, scale, 1000)
    time,I,S,R,times_to_sample,grid_I,grid_S=\
    Sellke_algo(tau,I_0,N,T_f,s_1,showplot=True,return_full=True)
    '''   
    
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    










