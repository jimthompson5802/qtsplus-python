

import numpy as np
from math import exp
from QtPyCommon import BasicQueue, QtsGlobalParameters,generate_time_vector, \
    BasicSimulator, side_by_side_barchart, ServerUnstable
from random import expovariate


class M_M_1(BasicQueue):
    """
    M/M/1 Model:  Poisson Arrivals to Single Exponential Server

    Poisson/exponential single-server queue with unlimited system capacity,
    FIFO:  calculates major measures of effectiveness and distributions
    for system size and wait time.
    GSTH, Section 2.2
    """

    def __init__(self,lam,st):
        """
        Parameters:
        
        lam: mean arrival rate
        st: mean service time
        """
        self.model_name="M/M/1"
        BasicQueue.__init__(self)

        self.lam = float(lam)
        self.st = float(st)
        self.c = 1
        self.K = np.Infinity

        
        #set up simulation parameters
        self.iat_fctn = expovariate
        self.st_fctn = expovariate

    def __setattr__(self,item,value):
        """
        Enforces data integrity of parameters between analytic
        and simulation parameters
        """
        #capture value
        self.__dict__[item] = value

        #enforce related values
        if item == "lam":
            self.__dict__["iat"] = 1/float(value)
            self.__dict__["iat_parms"] = [[float(value)]]
        elif item == "st":
            self.__dict__["mu"] = 1/float(value)
            self.__dict__["st_parms"] = [[1/float(value)]]

    def __print_model_specific_parameters__(self,type):
        #print model specific parameters
        if type == "analytic":
            self.print_parameters(["lam","iat","st","mu"])
        else:
            self.print_parameters(["iat_fctn","iat_parms",
                                  "st_fctn","st_parms"])

    def __print_model_specific_results__(self,type):
        #print model specific results
        if type == "analytic":
            self.an.print_results(["rho","W","Wq","L","Lq"])
        else:
            self.sim.print_results(["end_simulation_time",
                        "transactions_arrived",
                        "transactions_completed",
                        "avg_iat","var_iat","avg_st","var_st","rho",
                        "W","var_W","Wq","var_Wq",
                        "L","Lq"])

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
        self.an.traffic_intensity = self.lam / self.mu
        self.an.rho = self.an.traffic_intensity          #server utilization
        if self.an.rho >= 1.0:
            raise ServerUnstable("rho greater or equal to 1.0, rho=%g" % self.an.rho)
        self.an.L = self.an.rho/(1-self.an.rho)         #System size
        self.an.Lq = self.an.L - self.an.rho            #Queue length
        self.an.W = self.an.L/self.lam               #System wait time
        self.an.Wq = self.an.W - self.st             #Wait in queue

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
        rho_vec = [self.an.rho]*(max_n+1)
        self.an.Pn = np.column_stack([n_vec,map(Pn,n_vec,rho_vec),np.zeros(len(n_vec))])

        #calculate CDF or P[N=n]
        self.an.Pn[0,2] = self.an.Pn[0,1]
        for i in xrange(1,len(self.an.Pn)):
            self.an.Pn[i,2] = self.an.Pn[i,1] + self.an.Pn[i-1,2]

        
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
        t_vec = generate_time_vector(max_t,QtsGlobalParameters.max_plot_points)

        #define functions for W(t) and Wq(t)
        W_t = lambda t,lam,mu: 1.0 - exp(-(mu - lam)*t)
        Wq_t = lambda t,lam,mu: 1.0 - (lam/mu)*exp(-(mu - lam)*t)

        #create lam and mu vectors
        lam_vec = [self.lam] * len(t_vec)
        mu_vec = [self.mu] * len(t_vec)

        #compute CDF matrix
        self.an.CDFWt = np.column_stack([t_vec,map(W_t,t_vec,lam_vec,mu_vec)])
        self.an.CDFWqt = np.column_stack([t_vec,map(Wq_t,t_vec,lam_vec,mu_vec)])

    def simulate(self,maxNumber=1000,max_n=10,max_t=10):
        
        #set up simulation parameters for arrival and service times
        self.simulation_maxNumber = maxNumber

        #set up simulation object
        self.qsim =  BasicSimulator(self.sim,
                          iat_fctn = self.iat_fctn,
                          iat_parms = self.iat_parms,
                          st_fctn = self.st_fctn,
                          st_parms = self.st_parms)

        #run simulation model
        self.qsim.run_simulation(self.simulation_maxNumber,max_n,max_t)
        
if __name__ == "__main__":
    QtsGlobalParameters.simulation_checkpoint = 25000
    myq = M_M_1(3.0,0.2)
    myq.solve()
    myq.print_results("analytic")
    myq.plot_Pn(max_n=20,show_model_parameters=True,n_incr=5)
    myq.plot_CDF(show_model_parameters=True)
    myq.simulate(200000,20)
    myq.print_results("simulation")
    myq.plot_Pn("simulation",n_incr=5)
    myq.plot_CDF("simulation")
    width = 0.35
    x1 = myq.an.Pn[:,0]
    y1 = myq.an.Pn[:,1]
    x2 = myq.sim.Pn[:,0]
    y2 = myq.sim.Pn[:,1]


    side_by_side_barchart(x1,y1,x2,y2,label1="analytic",label2="simulation",n_incr=5)
    #plt.title('CDF Plot for analytic and simulation results')
    
    M_M_1(1,3).solve()
    
    
    
    



    
    
