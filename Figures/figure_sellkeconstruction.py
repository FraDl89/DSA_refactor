#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 15:28:39 2021

@author: fra
"""

import numpy as np
import matplotlib.pyplot as plt
import myPlotConfigs
import random
from scipy import stats
from matplotlib.text import OffsetFrom
import sys
sys.path.append('../')

from Selkealgo import Sellke_algo

'''
This produces the figure to explain how the sellke construction works
'''
if __name__=="__main__":
    N = 500
    beta =1.3
    seed = 12
    T_f = 10
    np.random.seed(seed)
    random.seed(seed)

    c = 1.9
    scale = 1                
    T_recover = stats.weibull_min.rvs(c,scale=scale,size=N)
    
    
    I_0 = 1
    tau = beta/(N-1)
    #epidemic parameters
    time,I,S,R,times_to_sample,grid_I,grid_S=\
    Sellke_algo(tau,I_0,N,T_f,T_recover,showplot=False,return_full=True)
    
    
    
    numb_points = 16
    cumulative_infectious_pressure = np.zeros(numb_points)
    
    
    for t in range(1,numb_points):
        cumulative_infectious_pressure[t] = cumulative_infectious_pressure[t-1]+I[t-1]*(time[t]-time[t-1])
        
        
    fig, ax1 = plt.subplots()    
    
    ax1.plot(time[:numb_points],cumulative_infectious_pressure, color='g')
    ax1.hlines(cumulative_infectious_pressure[-1],time[numb_points-1],4.4, colors='g')
    
    ax1.set_xlim(0,4.4)
    ax1.set_ylim(0)
    ax1.set_xlabel(r"Time",size=12)
    ax1.set_ylabel(r"Thresholds",size=12)
    #ax1.hlines(cumulative_infectious_pressure[1:numb_points-1],time[1:numb_points-1], time[2:numb_points])


    coordinate_base = (0.1,6)
    coordinate_base_shift=(0.3,6)
    p = -3
    ax1.annotate("Infections",fontsize=14,xy=coordinate_base,xytext=(0, 10), textcoords='offset points')

    epsilon = 0
    arrowprops=dict(facecolor='black',headwidth=4,width=1, shrink=0.007)
    ax1.annotate("",xy=(time[1],cumulative_infectious_pressure[1]+epsilon),
                 xycoords='data',xytext=coordinate_base_shift,textcoords='data',arrowprops=arrowprops)
    ax1.annotate("",xy=(time[3],cumulative_infectious_pressure[3]+epsilon),
                 xycoords='data',xytext=coordinate_base_shift,textcoords='data',arrowprops=arrowprops)
    ax1.annotate("",xy=(time[4],cumulative_infectious_pressure[4]+epsilon),
                 xycoords='data',xytext=coordinate_base_shift,textcoords='data',arrowprops=arrowprops)
    ax1.annotate("",xy=(time[5],cumulative_infectious_pressure[5]+epsilon),
                 xycoords='data',xytext=coordinate_base_shift,textcoords='data',arrowprops=arrowprops)
    ax1.annotate("",xy=(time[9],cumulative_infectious_pressure[9]+epsilon),
                 xycoords='data',xytext=coordinate_base_shift,textcoords='data',arrowprops=arrowprops)
    ax1.annotate("",xy=(time[11],cumulative_infectious_pressure[11]+epsilon),
                 xycoords='data',xytext=coordinate_base_shift,textcoords='data',arrowprops=arrowprops)
    ax1.annotate("",xy=(time[13],cumulative_infectious_pressure[13]),
                 xycoords='data',xytext=coordinate_base_shift,textcoords='data',arrowprops=arrowprops)


    ax1.hlines(cumulative_infectious_pressure[0],time[0],time[2], color='orange')
    ax1.vlines(time[2],cumulative_infectious_pressure[0],cumulative_infectious_pressure[2], color='k')

    ax1.hlines(cumulative_infectious_pressure[1],time[1],time[6], color='orange')
    ax1.vlines(time[6],cumulative_infectious_pressure[1],cumulative_infectious_pressure[6], color='k')

    ax1.hlines(cumulative_infectious_pressure[3],time[3],time[7], color='orange')
    ax1.vlines(time[7],cumulative_infectious_pressure[3],cumulative_infectious_pressure[7], color='k')


    ax1.hlines(cumulative_infectious_pressure[4],time[4],time[8], color='orange')
    ax1.vlines(time[8],cumulative_infectious_pressure[4],cumulative_infectious_pressure[8], color='k')


    ax1.hlines(cumulative_infectious_pressure[5],time[5],time[10], color='orange')
    ax1.vlines(time[10],cumulative_infectious_pressure[5],cumulative_infectious_pressure[10], color='k')

    ax1.hlines(cumulative_infectious_pressure[9],time[9],time[12], color='orange')
    ax1.vlines(time[12],cumulative_infectious_pressure[9],cumulative_infectious_pressure[12], color='k')
    
    ax1.hlines(cumulative_infectious_pressure[11],time[11],time[14], color='orange')
    ax1.vlines(time[14],cumulative_infectious_pressure[11],cumulative_infectious_pressure[14], color='k')
    
    ax1.hlines(cumulative_infectious_pressure[13],time[13],time[15], color='orange')
    ax1.vlines(time[15],cumulative_infectious_pressure[13],cumulative_infectious_pressure[15], color='k')
        
    left, bottom, width, height = [0.7, 0.17, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    
    space = np.linspace(0,5,100)
    ax2.plot(space,stats.weibull_min.pdf(space,c=c,scale=scale),color='orange')
    ax2.set_xticks([0,2,4])
    ax2.set_yticks([0,0.3,0.6])
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.set_title("Infectious period distribution", fontsize=9)
    import seaborn as sns
    sns.despine()
    
    
    ax1.spines['bottom'].set_color('#A3238E')
    ax1.spines['top'].set_color('#A3238E') 
    ax1.spines['right'].set_color('#A3238E')
    ax1.spines['left'].set_color('#A3238E')

    ax2.spines['bottom'].set_color('#A3238E')
    ax2.spines['top'].set_color('#A3238E') 
    ax2.spines['right'].set_color('#A3238E')
    ax2.spines['left'].set_color('#A3238E')

    plt.savefig("sellkeconstruction.eps",format='eps')
    '''
    plt.annotate("infection", xy=(time[1],cumulative_infectious_pressure[1]),
                 xycoords='data',xytext=(-70, 100), 
                 textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.05),)

    plt.annotate("infection", xy=(time[3],cumulative_infectious_pressure[3]),
                 xycoords='data',xytext=(-70, 100), 
                 textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.05),)



    plt.annotate("1st recovery", xy=(time[2],cumulative_infectious_pressure[2]/2),
                 xycoords='data',xytext=(40, -10), 
                 textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.05),)

    '''

