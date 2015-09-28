import numpy as np
from math import exp
import matplotlib.pyplot as plt
from QtPy.QtPyCommon import *
from random import expovariate


class M_M_1(BasicQueue):
    """
    M/M/1 Model:  Poisson Arrivals to Single Exponential Server

    Poisson/exponential single-server queue with unlimited system capacity,
    FIFO:  calculates major measures of effectiveness and distributions
    for system size and wait time.

    """

    def __init__(self,lam,st):
        """
        Parameters:
        
        lam: mean arrival rate
        st: mean service time
        """
        self.model_name="M/M/1"
        BasicQueue.__init__(self,lam,st)
 
    def __print_model_specific_parameters__(self):
        #print model specific parameters
        pass

    def __print_model_specific_results__(self):
        #print model specific results
        pass

    def solve(self,max_n=10,max_t=10):
        """
        Calculates basic performance results

        Returns:
        rho: server utilization
        W:  mean time at server
        Wq: mean wait time in queue
        L:  mean number of customers in system
        Lq: mean number of customers in queue

        """
        #calculate basic performance measures
        BasicQueue.update_basic_results(self)   #required...do not remove
        
        self.rho = self.traffic_intensity          #server utilization
        if self.rho >= 1:
            raise ServerUnstable("rho greater or equal to 1.0, rho=%g" % self.rho)
        self.L = self.rho/(1-self.rho)         #System size
        self.Lq = self.L - self.rho            #Queue length
        self.W = self.L/self.lam               #System wait time
        self.Wq = self.W - self.st             #Wait in queue

        self.calc_Pn(max_n)
        self.calc_CDF(max_t)

    def calc_Pn(self,max_n=10):
        """
        Calculates probability of N=n in the system

        Parameter:
        max_n: maximum system size for plotting

        Returns:
        Pn:  matrix of [n,P[N=n]]
        """
        Pn = lambda n,rho: (1.0 - rho)*rho**n
        n_vec = range(0,max_n+1)
        rho_vec = [self.rho]*(max_n+1)
        self.Pn = np.column_stack([n_vec,map(Pn,n_vec,rho_vec),np.zeros(len(n_vec))])

        #calculate CDF or P[N=n]
        self.Pn[0,2] = self.Pn[0,1]
        for i in xrange(1,len(self.Pn)):
            self.Pn[i,2] = self.Pn[i,1] + self.Pn[i-1,2]

        
    def calc_CDF(self,max_t=10):
        """
        Calculates CDF for W(t) and Wq(t)

        Parameter:
        max_t: maximum time for W(t) and Wq(t) plots

        Creates matrix:
        CDFWt: [t,W(t)]
        CDFWqt: [t,Wq(t)]
        """
        #create vector of time points (t)
        t_vec = generate_time_vector(max_t,GlobalParameters.max_plot_points)

        #define functions for W(t) and Wq(t)
        W_t = lambda t,lam,mu: 1.0 - exp(-(mu - lam)*t)
        Wq_t = lambda t,lam,mu: 1.0 - (lam/mu)*exp(-(mu - lam)*t)

        #create lam and mu vectors
        lam_vec = [self.lam] * len(t_vec)
        mu_vec = [self.mu] * len(t_vec)

        #compute CDF matrix
        self.CDFWt = np.column_stack([t_vec,map(W_t,t_vec,lam_vec,mu_vec)])
        self.CDFWqt = np.column_stack([t_vec,map(Wq_t,t_vec,lam_vec,mu_vec)])

    def simulate(self,maxNumber=1000,max_n=10,max_t=10):
        BasicQueue.update_basic_results(self)  #required...do not remove
        
        #set up simulation parameters for arrival and service times
        iat_fctn = expovariate
        iat_parms = [[self.lam]]
        st_fctn = expovariate
        st_parms = [[self.mu]]
        self.simulation_maxNumber = maxNumber

        #set up simulation object
        self.qsim =  BasicSimulator(iat_fctn = iat_fctn,
                          iat_parms = iat_parms,
                          st_fctn = st_fctn,
                          st_parms = st_parms)

        #run simulation model
        self.qsim.run_simulation(self.simulation_maxNumber,max_n,max_t)
        


if __name__ == "__main__":
    QtPyGlobalParameters.simulation_checkpoint = 25000
    myq = M_M_1(3.0,0.2)
    myq.solve()
    myq.print_results()
    #myq.plot_Pn(50,show_model_parameters=True,n_incr=5)
    #myq.plot_CDF(show_model_parameters=True)
    myq.simulate(200000,15)
    myq.plot_CDF();myq.plot_simulation_CDF()
    



    
    
