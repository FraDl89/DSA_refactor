#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 09:16:00 2020

@author: ld288
"""
import numpy as np
from scipy.interpolate import interp1d

from pyDOE import lhs

from PDE_solver import SIR_PDEroutine
from scipy.optimize import minimize
from scipy.signal import convolve

import pyswarms as ps



"""
A class that includes all the methods to compute the complete likelihood
using dynamic survival analysis, as described in the paper.

It is not necessary to input all the contributions from all the sources
of data. If some is missing, simply not inputting it will result in the
respective likelihood to be always 1. 
"""
class log_likelihood_models():
    '''
    Parameters
    ----------
    CIdist : function
        Infection hazard function
    ngrid: int
        grid resolution, needed only to solve the PDE, otherwise not necessary
    **kwargs: keywords arguments
        'T': float (default = 10.0)
            Final time, used to stop the PDE. Note that the PDE will be 
            solved on a grid np.linspace(0,T,ngrid).
        'hazard_rec': function
            Recovery hazard, used to compute ll_I
        'hazard_inf': function
            Infection hazard
        'rec_distr': function
            Recovery distribution, used to compute ll_R1
        'rho': float or "infer"
            Initial condition of the PdE. if "infer" it will be inferred from 
            data
        'rec_parms': numpy.array(dtype=float) or "infer"
            Recovery distribution parameters, if "infer" they will be inferred
            from data.
        'infect_times': numpy.array(dtype=float)
            times of infections, data used to compute ll_I
        'infectious_ages':numpy.array(dtype=float)
            infectious periods, data used to compute ll_R1
        'recov_times': numpy.array(dtype=float)
            recovery times, data used to compute ll_R2
        'hazard_inf_par': int (default = 1)
            number of parameters of the infectious distribution
            (Note, could be dropped in favour of something like len(inf_distr.__code__.co_varnames))
        'infer_recovery: boolean (default = True)
            whether infer recovery distribution
        'infer_infection: boolean (default = True)
            whether infer infectiousness distribution
            
    '''  
    def __init__(self, ngrid, **kwargs):        

        #Maybe move initial condition of the pde here?

        
        self.l_I=0
        self.l_R1 = 0
        self.l_R2 = 0
        self.nUgrid = ngrid
        self.nTgrid = self.nUgrid
        self.T_f = kwargs.get('T', 10)


        self.rho = kwargs.get("rho", "infer")
        
        self.hazard_rec = kwargs.get('hazard_rec',0)
        self.CIdist = kwargs.get('hazard_inf',0)
        self.rec_distr = kwargs.get('rec_distr',0)



        self.rec_parms = kwargs.get("rec_parms",1)
        self.CIdist_par = kwargs.get('hazard_inf_par',1)
      
        self.infect_times = kwargs.get('infect_times',0)
        self.recov_times = kwargs.get('recov_times',0)
        self.infectious_ages = kwargs.get('infectious_ages',0)

        self.infer_recovery  = kwargs.get('infer_recovery',True)
        self.infer_infection = kwargs.get('infer_infection',True)
        


    def pde(self, u_rec):
        #This is an auxiliary function that returns the pde object to be run in
        #The likelihood functions
        u=np.zeros(1+self.CIdist_par+self.rec_parms)
        if self.rho=='infer':
            u[0]=u_rec[0]
        else:
            u[0]=self.rho
        if self.infer_infection == True:
            for i in range(self.CIdist_par):
                if self.rho=='infer':
                    u[i+1]=u_rec[i+1]
                else:
                    u[i+1]=u_rec[i]
        else:
            u[1:self.CIdist_par+1]=self.CIdist_par
            
        if  self.infer_recovery == True:
            if self.rho=='infer':
                if  self.infer_infection == True:
                    u[1+self.CIdist_par:]=u_rec[1+self.CIdist_par:]
                else:
                    u[1+self.CIdist_par:]=u_rec[1:]
            else:
                if  self.infer_infection == True:
                    u[1+self.CIdist_par:]=u_rec[self.CIdist_par:]
                else:
                    u[1+self.CIdist_par:]=u_rec
 
        pdeObj_1 = SIR_PDEroutine(u[0], CIdist=self.CIdist, 
                                  CIdistParms=[*u[1:self.CIdist_par+1]], 
                                  recovDist=self.hazard_rec, 
                                  recovDistParms=u[self.CIdist_par+1:], 
                                  nTgrid=self.nTgrid, 
                                  nUgrid=self.nUgrid, T=self.T_f)
    
        self.initial_cond=np.exp(-pdeObj_1.tgrids)

        rec_times = self.rec_distr(pdeObj_1.tgrids,*u[self.CIdist_par+1:])                
    
                       
        return pdeObj_1,rec_times

    def ll_I(self,u_rec):
        '''
        log_likelihood to infer recovery distribution parameters from infection times
        This is l_I^1 in the paper
        
        Parameters:
            urec -> vector, parameters of the recovery hazard to infer
        Returns: log likelihood of l_I. Note, if no data given, will return 0
        '''
        if type(self.infect_times)==int:
            return 0
        else:

            
            xS = interp1d(self.pdeObj_1.tgrids, self.X)
            YtsI= interp1d(self.pdeObj_1.tgrids, self.Yts)
            
            infect_times = self.infect_times[self.infect_times<self.T_f]

            S_eval = np.log(xS(infect_times)) 
            Yts_eval = np.log(YtsI(infect_times)) 
            tau_T = np.log(1 - xS(infect_times[-1]))
            
            #print((np.log(u_rec[1])- tau_T)*len(infect_times) + np.sum( S_eval + I_eval))
            return ((- tau_T)*len(infect_times) + np.sum( S_eval + Yts_eval))

        
    def ll_R1(self,u_rec):
        '''
        log_likelihood to infer recovery distribution parameters from infectious ages
        This is ll_R1 in the paper
        
        Parameters:
            urec -> vector, parameters of the recovery distribution to infer

        Returns: log likelihood of l_R1. Note, if no data given, will return 0
        '''
        if type(self.infectious_ages)==int:
            return 0
        else:
            u=np.zeros(1+self.CIdist_par+self.rec_parms)
            if self.rho=='infer':
                u[0]=u_rec[0]
            else:
                u[0]=self.rho
            if self.infer_infection == True:
                for i in range(self.CIdist_par):
                    if self.rho=='infer':
                        u[i+1]=u_rec[i+1]
                    else:
                        u[i+1]=u_rec[i]
            else:
                u[1:self.CIdist_par+1]=self.CIdist_par
                
            if  self.infer_recovery == True:
                if self.rho=='infer':
                    if  self.infer_infection == True:
                        u[1+self.CIdist_par:]=u_rec[1+self.CIdist_par:]
                    else:
                        u[1+self.CIdist_par:]=u_rec[1:]
                else:
                    if  self.infer_infection == True:
                        u[1+self.CIdist_par:]=u_rec[self.CIdist_par:]
                    else:
                        u[1+self.CIdist_par:]=u_rec

             
        return (np.sum(np.log(self.rec_distr(self.infectious_ages,u[self.CIdist_par+1:]))))





    def ll_R2(self,u_rec):
        '''
        log_likelihood to infer recovery distribution parameters from recovery times
        This is ll_R2 in the paper
        
        Parameters:
            urec -> vector, parameters of the recovery hazard to infer

        Returns: log likelihood of l_R2. Note, if no data given, will return 0
        '''
        #To do, find a way to compute the convolution-form.
        if type(self.recov_times)==int:
            return 0
        else:

            
            recov_times = self.recov_times[self.recov_times<self.T_f]



            
 
            ft = self.X*self.Yts
            
            gt = convolve(ft,self.rec_times)[1:len(self.pdeObj_1.tgrids)]*self.dx
            gt = gt/(np.sum(gt)*self.dx)
            #print(min(ft),max(ft), min(gt), max(gt))
            loggt = np.log(gt)
            res=interp1d(self.pdeObj_1.tgrids[1:], loggt)(recov_times[recov_times>self.pdeObj_1.tgrids[1]]).sum()
            return res

    def evaluate_likelihood(self,x):
            self.pdeObj_1,self.rec_times = self.pde(x)
            self.X,self.Y, self.Yts=self.pdeObj_1.finDiffUpdate(Yts=True,initial_cond=self.initial_cond)
            self.dx = self.pdeObj_1.dx
                       
            return -(self.ll_I(x) + self.ll_R1(x) +self.ll_R2(x))    



    def minimize_likelihood(self,lb,ub,use_pyswarm=True, maxiter=100,swarmsize=100, c1=1.7, c2=1.7, w=0.6):
        '''
        Parameters
        ----------
        lb : array like,float
            lower bounds array. Len(lb) is the number of parameters to infer
        ub : array like,float
            upper bounds array. Len(ub) is the number of parameters to infer
        initial_points: INT
            number of random samples of the region delimited by the bounds to start 
            the maximiser with
        Returns
        -------
        float
            result of minimization.

        '''
  

        def neg_log_likelihood(u):

            results = np.zeros(len(u))
            for index,x in enumerate(u):
  
                self.pdeObj_1,self.rec_times = self.pde(x)
                self.X,self.Y, self.Yts=self.pdeObj_1.finDiffUpdate(Yts=True,initial_cond=self.initial_cond)
                self.dx = self.pdeObj_1.dx
                results[index]=-(self.ll_I(x) + self.ll_R1(x) +self.ll_R2(x))
          
            return results
            
        if use_pyswarm==True:
            options = {'c1':c1, 'c2': c2, 'w':w}
            bounds = np.array([(lb[i],ub[i]) for i in range(len(lb))])
   
            optimizer = ps.single.GlobalBestPSO(n_particles=swarmsize, dimensions=len(lb), options=options,bounds=[lb,ub])

            fopt, xopt = optimizer.optimize(neg_log_likelihood, iters=maxiter)

            #fopt, xopt = pso(neg_log_likelihood,lb,ub,maxiter=maxiter,swarmsize=swarmsize,minfunc=1e-4, debug=False)
            
            #bounds = np.array([(lb[i],ub[i]) for i in range(len(lb))])
            #res=minimize(neg_log_likelihood, xopt, bounds = bounds, options={'maxiter':550,},method='TNC')
            #fun = shgo(neg_log_likelihood, bounds, options={'ftol':1e-5})
            return xopt, fopt
        else:
            bounds = np.array([(lb[i],ub[i]) for i in range(len(lb))])
            #Latin hypercube sampling to generate a uniform distribution of points
              
            
            #The idea is to launch the minimizer on a grid and find the value where the
            #negative log likelihood is lower. This is done to avoid incurring in local
            #minima.    
            lhd = lhs(len(lb), samples=5000)*(ub-lb)+lb    
            
            evaluation = np.PINF
            x = 0
            for x0 in lhd:
                result = neg_log_likelihood(x0)
                if result < evaluation and np.isfinite(result):
                    x = x0
                    evaluation = result
            #print(evaluation,x)
            res = minimize(neg_log_likelihood, x, bounds = bounds, options={'maxiter':200,},method='TNC')

        return res
        

       

"""
A simpler class that includes all the methods to compute the complete likelihood
using dynamic survival analysis, as described in the paper. This is valid 
ONLY if the infectiousness is exponential distributed (as for synthetic data)

It is not necessary to input all the contributions from all the sources
of data. If some is missing, simply not inputting it will result in the
respective likelihood to be always 1. 
"""
class log_likelihood_models_simple():
    '''
    Parameters
    ----------
    CIdist : function
        Infection hazard function
    ngrid: int
        grid resolution, needed only to solve the PDE, otherwise not necessary
    **kwargs: keywords arguments
        'T': float (default = 10.0)
            Final time, used to stop the PDE. Note that the PDE will be 
            solved on a grid np.linspace(0,T,ngrid).
        'hazard_rec': function
            Recovery hazard, used to compute ll_I
        'hazard_inf': function
            Infection hazard
        'rec_distr': function
            Recovery distribution, used to compute ll_R1
        'rho': float or "infer"
            Initial condition of the PdE. if "infer" it will be inferred from 
            data
        'rec_parms': numpy.array(dtype=float) or "infer"
            Recovery distribution parameters, if "infer" they will be inferred
            from data.
        'infect_times': numpy.array(dtype=float)
            times of infections, data used to compute ll_I
        'infectious_ages':numpy.array(dtype=float)
            infectious periods, data used to compute ll_R1
        'recov_times': numpy.array(dtype=float)
            recovery times, data used to compute ll_R2
        'hazard_inf_par': int (default = 1)
            number of parameters of the infectious distribution
            (Note, could be dropped in favour of something like len(inf_distr.__code__.co_varnames))
        'infer_recovery: boolean (default = True)
            whether infer recovery distribution
        'infer_infection: boolean (default = True)
            whether infer infectiousness distribution
            
    '''  
    def __init__(self, ngrid, **kwargs):        

        #Maybe move initial condition of the pde here?

        
        self.l_I=0
        self.l_R1 = 0
        self.l_R2 = 0
        self.nUgrid = ngrid
        self.nTgrid = self.nUgrid
        self.T_f = kwargs.get('T', 10)


        self.rho = kwargs.get("rho", "infer")
        
        self.hazard_rec = kwargs.get('hazard_rec',0)
        self.CIdist = kwargs.get('hazard_inf',0)
        self.rec_distr = kwargs.get('rec_distr',0)



        self.rec_parms = kwargs.get("rec_parms",1)
        self.CIdist_par = kwargs.get('hazard_inf_par',1)
      
        self.infect_times = kwargs.get('infect_times',0)
        self.recov_times = kwargs.get('recov_times',0)
        self.infectious_ages = kwargs.get('infectious_ages',0)

        self.infer_recovery  = kwargs.get('infer_recovery',True)
        self.infer_infection = kwargs.get('infer_infection',True)
        


    def pde(self, u_rec):
        #This is an auxiliary function that returns the pde object to be run in
        #The likelihood functions
        u=np.zeros(1+self.CIdist_par+self.rec_parms)
        if self.rho=='infer':
            u[0]=u_rec[0]
        else:
            u[0]=self.rho
        if self.infer_infection == True:
            for i in range(self.CIdist_par):
                if self.rho=='infer':
                    u[i+1]=u_rec[i+1]
                else:
                    u[i+1]=u_rec[i]
        else:
            u[1:self.CIdist_par+1]=self.CIdist_par
            
        if  self.infer_recovery == True:
            if self.rho=='infer':
                if  self.infer_infection == True:
                    u[1+self.CIdist_par:]=u_rec[1+self.CIdist_par:]
                else:
                    u[1+self.CIdist_par:]=u_rec[1:]
            else:
                if  self.infer_infection == True:
                    u[1+self.CIdist_par:]=u_rec[self.CIdist_par:]
                else:
                    u[1+self.CIdist_par:]=u_rec
 
        pdeObj_1 = SIR_PDEroutine(u[0], CIdist=self.CIdist, 
                                  CIdistParms=[*u[1:self.CIdist_par+1]], 
                                  recovDist=self.hazard_rec, 
                                  recovDistParms=u[self.CIdist_par+1:], 
                                  nTgrid=self.nTgrid, 
                                  nUgrid=self.nUgrid, T=self.T_f)
    
        self.initial_cond=np.exp(-pdeObj_1.tgrids)

        rec_times = self.rec_distr(pdeObj_1.tgrids,*u[self.CIdist_par+1:])                
    
                       
        return pdeObj_1,rec_times

    def ll_I(self,u_rec):
        '''
        log_likelihood to infer recovery distribution parameters from infection times
        This is l_I^1 in the paper
        
        Parameters:
            urec -> vector, parameters of the recovery hazard to infer
        Returns: log likelihood of l_I. Note, if no data given, will return 0
        '''
        if type(self.infect_times)==int:
            return 0
        else:

            
            xS = interp1d(self.pdeObj_1.tgrids, self.X)
            YtsI= interp1d(self.pdeObj_1.tgrids, self.Yts)
            
            infect_times = self.infect_times[self.infect_times<self.T_f]

            S_eval = np.log(xS(infect_times)) 
            Yts_eval = np.log(YtsI(infect_times)) 
            tau_T = np.log(1 - xS(infect_times[-1]))
            
            #print((np.log(u_rec[1])- tau_T)*len(infect_times) + np.sum( S_eval + I_eval))
            return ((- tau_T)*len(infect_times) + np.sum( S_eval + Yts_eval))

        
    def ll_R1(self,u_rec):
        '''
        log_likelihood to infer recovery distribution parameters from infectious ages
        This is ll_R1 in the paper
        
        Parameters:
            urec -> vector, parameters of the recovery distribution to infer

        Returns: log likelihood of l_R1. Note, if no data given, will return 0
        '''
        if type(self.infectious_ages)==int:
            return 0
        else:
            u=np.zeros(1+self.CIdist_par+self.rec_parms)
            if self.rho=='infer':
                u[0]=u_rec[0]
            else:
                u[0]=self.rho
            if self.infer_infection == True:
                for i in range(self.CIdist_par):
                    if self.rho=='infer':
                        u[i+1]=u_rec[i+1]
                    else:
                        u[i+1]=u_rec[i]
            else:
                u[1:self.CIdist_par+1]=self.CIdist_par
                
            if  self.infer_recovery == True:
                if self.rho=='infer':
                    if  self.infer_infection == True:
                        u[1+self.CIdist_par:]=u_rec[1+self.CIdist_par:]
                    else:
                        u[1+self.CIdist_par:]=u_rec[1:]
                else:
                    if  self.infer_infection == True:
                        u[1+self.CIdist_par:]=u_rec[self.CIdist_par:]
                    else:
                        u[1+self.CIdist_par:]=u_rec

             
        return (np.sum(np.log(self.rec_distr(self.infectious_ages,u[self.CIdist_par+1:]))))





    def ll_R2(self,u_rec):
        '''
        log_likelihood to infer recovery distribution parameters from recovery times
        This is ll_R2 in the paper
        
        Parameters:
            urec -> vector, parameters of the recovery hazard to infer

        Returns: log likelihood of l_R2. Note, if no data given, will return 0
        '''
        #To do, find a way to compute the convolution-form.
        if type(self.recov_times)==int:
            return 0
        else:

            
            recov_times = self.recov_times[self.recov_times<self.T_f]



            
 
            ft = self.X*self.Yts
            
            gt = convolve(ft,self.rec_times)[1:len(self.pdeObj_1.tgrids)]*self.dx
            gt = gt/(np.sum(gt)*self.dx)
            
            loggt = np.log(gt)
            res=interp1d(self.pdeObj_1.tgrids[1:], loggt)(recov_times[recov_times>self.pdeObj_1.tgrids[1]]).sum()
            return res



    def evaluate_likelihood(self,x):
            self.pdeObj_1,self.rec_times = self.pde(x)
            self.X,self.Y, self.Yts=self.pdeObj_1.finDiffUpdate(Yts=True,initial_cond=self.initial_cond)
            self.dx = self.pdeObj_1.dx
                       
            return -(self.ll_I(x) + self.ll_R1(x) +self.ll_R2(x))        




    def minimize_likelihood(self,lb,ub,use_pyswarm=True, maxiter=100,swarmsize=100):
        
        from pyswarm import pso
        #This is pyswarm, a very simple particle swarm algorithm whose implementation
        #can be found at https://pypi.org/project/pyswarm/
        #It is simpler to use than pyswarms, but does the trick in this case
        '''

        Parameters
        ----------
        lb : array like,float
            lower bounds array. Len(lb) is the number of parameters to infer
        ub : array like,float
            upper bounds array. Len(ub) is the number of parameters to infer
        initial_points: INT
            number of random samples of the region delimited by the bounds to start 
            the maximiser with
        Returns
        -------
        float
            result of minimization.

        '''
  

        def neg_log_likelihood(x):
            self.pdeObj_1,self.rec_times = self.pde(x)
            self.X,self.Y, self.Yts=self.pdeObj_1.finDiffUpdate(Yts=True,initial_cond=self.initial_cond)
            self.dx = self.pdeObj_1.dx
                       
            return -(self.ll_I(x) + self.ll_R1(x) +self.ll_R2(x))

        if use_pyswarm==True:
               xopt,fopt = pso(neg_log_likelihood,lb,ub,maxiter=maxiter,swarmsize=swarmsize,minfunc=1e-4, debug=False)
               bounds = np.array([(lb[i],ub[i]) for i in range(len(lb))])
               print(xopt,fopt) 
               res=minimize(neg_log_likelihood, xopt, bounds = bounds, options={'maxiter':550,},method='TNC')
               #fun = shgo(neg_log_likelihood, bounds, options={'ftol':1e-5})
               return res
        else:
            bounds = np.array([(lb[i],ub[i]) for i in range(len(lb))])
            #Latin hypercube sampling to generate a uniform distribution of points
              
            
            #The idea is to launch the minimizer on a grid and find the value where the
            #negative log likelihood is lower. This is done to avoid incurring in local
            #minima.    
            lhd = lhs(len(lb), samples=5000)*(ub-lb)+lb    
            
            evaluation = np.PINF
            x = 0
            for x0 in lhd:
                result = neg_log_likelihood(x0)
                if result < evaluation and np.isfinite(result):
                    x = x0
                    evaluation = result
            #print(evaluation,x)
            res = minimize(neg_log_likelihood, x, bounds = bounds, options={'maxiter':200,},method='TNC')

        return res



