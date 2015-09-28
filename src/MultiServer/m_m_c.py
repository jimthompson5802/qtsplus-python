from QtPy.QtPyCommon import *
#Other required imports
import numpy as np
from math import exp
from random import expovariate

class M_M_c(GlobalParameters,BasicQueue):
    """
    M/M/c Model:  M/M/c: Poisson Arrivals to Multiple Exponential Servers

    Poisson/exponential multi-server queue with unlimited system capacity,
    FIFO: calculates the major measures of effectiveness.
    GSTH, Section 2.3
    """

    def __init__(self,lam,st,c):
        """
        Parameters:
        -----------
        lam: mean arrival rate
        st: mean service time
        c: number of servers
        """
        
        self.model_name="M/M/c"
        BasicQueue.__init__(self)

        self.lam = float(lam)
        self.st = float(st)
        self.c = c
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
            self.print_parameters(["lam","iat","st","mu","c"])
        else:
            self.print_parameters(["iat_fctn","iat_parms",
                                  "st_fctn","st_parms","c"])

       

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


    def solve(self):
        """
        Calculates basic performance results

        Returns:
        --------
        rho: server utilization
        W:  mean time at server
        Wq: mean wait time in queue
        L:  mean number of customers in system
        Lq: mean number of customers in queue
        P0: Probability the queueing station is empty
        1-Wq(0): Prob. arriving customer is delayed in queue
        """
        
        #compute core model results
        self.an.traffic_intensity  = self.lam / self.mu
        self.an.rho = self.an.traffic_intensity / float(self.c)
        if self.an.rho >= 1:
            raise ServerUnstable("rho greater or equal to 1.0, rho=%g" % self.an.rho)
        self.an.P0 = M_M_c_P0(self.an.traffic_intensity,self.c)
        self.an.Lq = M_M_c_Lq(self.an.traffic_intensity,self.c)
        self.an.L = self.an.Lq + self.an.traffic_intensity
        self.an.W = self.an.L / self.lam
        self.an.Wq = self.an.W - self.st

    def calc_Pn(self,max_n):
        """
        Calculates probability of N=n in the system

        Parameter:
        max_n: maximum sys.tem size for plotting

        Returns:
        Pn:  matrix of [n,P[N=n]]
        """
        #create parameter vectors
        n_vec = range(0,max_n+1)
        Pn_vec = [0.0]*(max_n+1)

        #oompute P[N=n]
        for n in range(0,max_n+1):
            if n == 0:
                Pn = self.an.P0
            elif n < self.c:
                Pn *= self.an.traffic_intensity / n
            else:
                Pn *= self.an.traffic_intensity / self.c
            Pn_vec[n] = Pn

        self.an.Pn = np.column_stack([n_vec,Pn_vec])

    def calc_CDF(self,max_t):
        """
        Calculates CDF for Wq(t)

        Parameter:
        max_t: maximum time for Wq(t) plot

        Returns:
        CDF: matrix of [t,Wq(t)]
        """
        #create vector of time points (t)
        max_pts = GlobalParameters.max_plot_points
        incr = float(max_t)/max_pts
        n_vec = range(0,max_pts+1)
        t_vec = [incr] * (max_pts+1)
        t_vec = map(lambda n,t: n*t,n_vec,t_vec)

        #define functions for Wq(t)
        Wq_t = lambda t,lam,mu,c,Lq: \
            1.0 - (Lq * (1 - (lam/(c*mu))) / (lam/(c*mu))) * exp(-(c * mu - lam) * t)

        #create required paramenter vectors
        lam_vec = [self.lam] * (max_pts+1)
        mu_vec = [self.mu] * (max_pts+1)
        c_vec = [self.c] * (max_pts+1)
        Lq_vec = [self.an.Lq] * (max_pts+1)

        #compute CDF matrix
        self.an.CDFWqt = np.column_stack([t_vec,map(Wq_t,t_vec,lam_vec,mu_vec,c_vec,Lq_vec)])


    def simulate(self,max_simulate=1000,max_n=10,max_t=10):
        self.simulation_maxNumber = max_simulate

        #set up simulation object
        self.qsim =  BasicSimulator(self.sim,
                          iat_fctn = self.iat_fctn,
                          iat_parms = self.iat_parms,
                          st_fctn = self.st_fctn,
                          st_parms = self.st_parms,
                          c=self.c)

        #run simulation model
        self.qsim.run_simulation(self.simulation_maxNumber,max_n,max_t)
        

    
def M_M_c_P0(r,c):
    """
    Computes P[N=0] for M/M/c
    r: traffic intensity (lam/mu)
    c: number of servers
    """
    T1 = 1.0
    T2 = 1.0
    for n in range(1,c):
        T2 = T2 * (r / n)
        T1 = T1 + T2
    T2 = T2 * (r / c)

    return 1.0 / (T1 + T2 / (1.0 - (r / c)))

def M_M_c_Lq(r,c):
    """
    Computes Lq
    r: traffic intensity (lam/mu)
    c: number of servers
    """
    rho = r / c

    T1 = 1.0
    for n in range(1,c+1):
        T1 = T1 * r/n

    return ((T1 * rho) / (1.0 - rho) ** 2) * M_M_c_P0(r, c)


if __name__ == "__main__":
    QtPyGlobalParameters.simulation_checkpoint = 50000
    q = M_M_c(5,1,8)
    q.solve()
    q.print_results()
    q.simulate(200000,15,20)
    q.print_results("simulation")
    q.calc_Pn(15)
    q.plot_Pn("analytic",15,n_incr=5)
    q.plot_Pn("simulation",n_incr=5)
    q.calc_CDF(20)
    q.plot_CDF(show_model_parameters=True)
    q.plot_CDF("simulation")

    
    
