'''
Copyright 2020-2022 Francesco Di Lauro. All Rights Reserved.
See LICENCE file for details.
'''  

import numpy as np
import scipy.stats as stats

from scipy.sparse import diags


class SIR_PDEroutine():
    '''
    Class to solve the PDE delined in Survival Analysis Dynamics paper, of the
    form
    $$
    (\partial_t + \partial_s) y_I(t,s) = - \gamma(s) y_I(t,s)
    d/dt y_S(t) = - \beta y_S(t) \int_0^\infty y(t,u) du
    $$
    with the initial conditions such that
    $y_S(0) = 1$
    $\int_0^\infty y_I(0,s) ds = \rho$ 
    and boundary condition
    $y_I(t,0) =  \int_0^\infty y_S(t,s)*beta*I ds$
    
    Params:
        rho: float. number of infected/number of susceptible at time 0
        CIdist: function.  Hazard function of infectivity
        CIdistParms: vector of floats. Parameters of CIdist function
        recovDist: function.  Hazard function of recovery
        recovDistParms: vector of floats. Parameters of recovDist function
    kwargs:
        T: float, cut-off time
        nTgrid: float, number of intervals that divide (0,T), time-keeping
        nUgrid: float, number of intervals that divide (0,T), ages-keeping. 
        N.b. nUgrid = nTgrid at the moment.
    '''
    def __init__(self, rho, CIdist, CIdistParms, recovDist, recovDistParms, **kwargs):
        self.rho = rho
        self.CIdist = CIdist
        self.CIdistParms = CIdistParms
        self.recovDist= recovDist
        self.recovDistParms = recovDistParms
        if kwargs.get('T') == None:
            self.T = 10.0
            #set default cut-off time if not provided
        else:
            self.T = kwargs.get('T')
        if kwargs.get ('nTgrid') == None:
            self.nTgrid = 10**3
            #set default number of grids along t-axis if not provided
        else:
            self.nTgrid = kwargs.get('nTgrid')
        if kwargs.get('nUgrid') == None:
            self.nUgrid = 10**3
            #set default number of grids along the u-axis if not provided
        else:
            self.nUgrid = kwargs.get('nUgrid')
        u0 = 0.0
        t0 = 0.0
        #create the grids along the t- and the u- axes
        self.ugrids = np.linspace (u0, self.T, self.nUgrid)
        self.tgrids = np.linspace (t0, self.T, self.nTgrid)
        self.dx = np.diff(self.ugrids)[0]
        self.mp =  np.linspace(0, self.T, self.nUgrid-1) + self.dx / 2.0  # element midpoints where γ is evaluated

    def infection_hazard(self,u):
        #defines the hazard function for the contact intervals
        #the function beta in the manuscript

        return self.CIdist(u,*self.CIdistParms)
    def  recv_hazard(self,u):
        #defines the hazard function for the recovery distribution
        #the function gamma in the manuscript

        return self.recovDist(u,*self.recovDistParms)

    def totInfecPressure(self,yvals):
        #this is the sum of infection pressure over all infectious ages
        #sum over beta times proportion of infected individuals
        dh = self.ugrids[1] - self.ugrids[0]
        infec_pressure = np.sum ([self.infection_hazard (self.ugrids[iter]) * yvals[iter] * dh for iter in range (self.nUgrid)])
        if infec_pressure<0:
            return 0
        else:
            return infec_pressure

    def suscEstimate(self,i,y):
        dt = self.tgrids[1] - self.tgrids[0]
        temp = np.sum([dt * self.totInfecPressure(yvals=y[j]) for j in range(i+1)])
        # for j in range(i+1):
        #     yvals = y[j]
        #     temp += dt * self.totInfecPressure(yvals)
        susceptible = np.exp(-temp)
        return susceptible



    def finDiffUpdate(self, **kwargs):
        '''
         Solves the PDE, returns X[t] and I_(t,s)}
         kwargs: 
             Yts: boolean (default=False)
                 if True, returns also \int <\beta, Y(t,s)> ds
             incond: vector of length nTgrid, default 'Unif'
                 y(0,s), initial condition. If "Unif", it will be uniform(0,5) days
                 
        '''
        return_Yts = kwargs.get("Yts", "False")
        initialcond=kwargs.get("incond", "Unif")
        X = np.zeros([self.nTgrid])  # fraction of susceptible individuals

        Y = np.zeros_like(X)  # for the initial condition
        X_I = np.zeros([self.nTgrid])

        if initialcond=='Unif':
  
            #Uniform for 5 days
            oneday=self.nUgrid/self.T
            fiveday = int(5*oneday)
            
            Y[:fiveday] = self.rho / (fiveday*self.dx)      
            Y[0] = np.sum(self.infection_hazard(self.tgrids[:fiveday])*Y[:fiveday]*self.dx)
            if return_Yts==True:
                Yts=np.zeros([self.nTgrid])
                Yts[0] = Y[0]
        else:
            Y = self.rho*initialcond  #Order 0
            Y[0] = np.sum(self.infection_hazard(self.tgrids)*Y*self.dx)
            
            #Note, boundary condition is not necessarily satisfied at t=0
            #This could be implemented by iterative methods, but it does not
            #change noticeably the outcome of the simulation
            
            #nonzero=np.where(Y>=Y[0]) #First order
            #Y[nonzero] = Y[nonzero]-Y[0]/(len(Y[nonzero]))
            
            if return_Yts==True:
                Yts=np.zeros([self.nTgrid])            
                Yts=Y[0]
        
        X_I[0] = self.rho  #Initial number of susceptibles
        
        #These (A,A1,A2) are different approximation of the pde operator
        #A is exact (but it can be costly to compute)
        #A1 and A2 are first order approximations
        #expgamma = np.exp(-self.dx*self.recv_hazard(self.mp))
        # ss = np.linspace(0, S, nS)  # nodal positions
        #expgamma = np.exp(-self.dx*self.recv_hazard(self.mp))

        inf_haz_dx=self.infection_hazard(self.tgrids)* self.dx
        # Sparse CSR matrix to approximate PDE operator: explicit Semi-Lagrangian method
        #A1 = diags([1.0 - dx * gamma], [-1]).tocsr()
        #A2 = diags([1.0 / (1.0 + dx * gamma)], [-1]).tocsr()
        #A = diags([expgamma], [-1]).tocsr()
        A = diags([1.0 / (1.0 + self.dx * self.recv_hazard(self.mp))], [-1]).tocsr()
        #A = diags([expgamma], [-1]).tocsr()
        #A= diags([1.0 - self.dx * self.recv_hazard(self.mp)], [-1]).tocsr()
        X[0] = 1.0  # initial condition of x_S
        # Solve initial boundary value problems
        for t in range(self.nTgrid - 1):
            Y = A.dot(Y)  # PDE propagation
            intY = np.sum(Y*inf_haz_dx)  # intY[t] = ∫_0^∞ y(t,s) ds
            X[t + 1] = X[t] / (1.0 +  intY* self.dx)  # update X
            X_I[t+1] = np.sum(Y) * self.dx
            if return_Yts==True:
                Yts[t+1]=intY
            Y[0] = X[t + 1] * intY  # update Y at boundary with implicit scheme
        if return_Yts==True:
            return X, X_I, Yts
        else:
            return X,X_I




if __name__=="__main__":
        #If you run this bit of code, you get a couple of examples of how the
        #Pde works
    
        grids = 10000
 
        #Example: classical SIR (exponential infectious period and infectiousness)        

        def rec_haz(u, *recovDistParams):    
            beta = float(recovDistParams[0])
            return beta*np.ones_like(u)   

        def inf_distr(u,*CIdistParms):
            beta = float(CIdistParms[0])
            return beta*np.ones_like(u)   

            
        import matplotlib.pyplot as plt
        


        pde= SIR_PDEroutine(0.01, CIdist=inf_distr, CIdistParms=[2],\
                                  recovDist=rec_haz, recovDistParms=[1],\
                                      nTgrid=grids, nUgrid=grids, T=25)

        X,Y=pde.finDiffUpdate() #Initial condition here is uniform (not specified)

        Yi = np.diff(pde.tgrids)[0] * np.array([np.sum(Y[t]) for t in range(pde.nTgrid)])

        plt.plot(pde.tgrids,Y, color='red', label='Infected')
        plt.plot(pde.tgrids,X, color='green', label='Susceptible')
        plt.title('beta = exp(2), gamma = exp(1)')
        plt.xlabel('time')
        plt.legend()



        #Example: Gamma for recovery, weibull for infectiousness
            

        def rec_haz2(u, *recovDistParams):   
            #Gamma distribution in terms of mean and variance 
            
            a = float(recovDistParams[0])**2/float(recovDistParams[1])
            scale = float(recovDistParams[1])/float(recovDistParams[0])   
            tol = 1e-10
            #Basically: use de l'hopital when the ratio becomes 0/0
            #Otherwise go with definition. This regularises a lot the numerics
            
            x = np.where(stats.gamma.cdf(u,a=a,scale=scale)>1-tol,
                         1/scale - (a-1)/u,
                         stats.gamma.pdf(u,a=a,scale=scale)/(1- stats.gamma.cdf(u,a=a,scale=scale)))
            return x               
                
        def inf_distr2(u,*CIdistParms):
           #Weibull distribution hazard by definition
           shape=CIdistParms[0]
           scale=CIdistParms[1]
           return shape/scale * (u/scale)**(shape-1)


        

        pde2= SIR_PDEroutine(0.01, CIdist=inf_distr2, CIdistParms=[5,1.8],\
                                  recovDist=rec_haz2, recovDistParms=[1.5,0.5],\
                                      nTgrid=grids, nUgrid=grids, T=20)



            
        initialcondition=np.exp(-pde.tgrids)
        #Initial condition: infected individual time from infection is
        #exp(-s)

        X_2,Y_2=pde2.finDiffUpdate(initialcond=initialcondition)

        Yi_2 = np.diff(pde2.tgrids)[0] * np.array([np.sum(Y_2[t]) for t in range(pde2.nTgrid)])

        
        #Plot the infectiousness and recovery distribution
        plt.figure()

        x=np.linspace(0,10,400)
        infdist=stats.weibull_min(5,scale=1.8)
        
        a = 1.5**2/1.8
        scale= 1.8/1.5
        recdist=stats.gamma(a=a, scale=scale)
            
        plt.plot(x,recdist.pdf(x), color='orange', label='recovery distribution')        
        plt.plot(x,infdist.pdf(x), color='blue', label='infectiousness distribution')
        plt.xlabel('time')
        plt.ylabel('density')
        plt.title('Weib(5,1.8), gamma(1,25,1.2))')
        plt.legend()
        
        #Plot epidemic curves
        plt.figure()
        plt.plot(pde2.tgrids,Y_2, color='red', label='Infected')
        plt.plot(pde2.tgrids,X_2, color='green', label='Susceptible')
        plt.xlabel('time')
        plt.ylabel('solution')
        plt.title('Weib(5,1.8), gamma(1,25,1.2)')
        plt.legend()

