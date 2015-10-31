import numpy as np
from math import exp
from QtPyCommon import BasicQueue, QtsGlobalParameters,generate_time_vector, \
    BasicSimulator, side_by_side_barchart
from random import expovariate

class M_M_1_K(QtsGlobalParameters,BasicQueue):
    """
    M/M/1/K Model: Poisson Arrivals to a Space-Limited Single Exponential Server

    Poisson/exponential single-server queue with limit on system capaicty
    of K:  calculates major measures of effectiveness and distributions for
    system size and wait time.
    GSTH, Section 2.5
    """
    def __init__(self,lam,st,K):
        """
        Parameters;
        lam: mean arrival rate
        st: mean service time
        K: size of queue
        """
        self.model_name = "M/M/1/K"
        BasicQueue.__init__(self)

        self.lam = float(lam)
        self.st = float(st)
        self.c = 1
        self.K = K

        #set up simulation parameters
        self.iat_fctn = expovariate
        self.st_fctn = expovariate

        #register model specific results
        for obj in [self.an,self.sim]:
            obj.add_result_variable("TA","Expected number turned away/unit time(TA)=%g")
            obj.add_result_variable("lam_eff","Effective arrival rate(lam_eff)=%g")
            obj.add_result_variable("rho_eff","Effective server utilization(rho_eff)=%g")

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
            self.an.print_results(["rho","rho_eff","lam_eff","W","Wq","L","Lq"])
        else:
            self.sim.print_results(["end_simulation_time",
                        "transactions_arrived",
                        "transactions_completed",
                        "avg_iat","var_iat","avg_st","var_st","rho",
                        "W","var_W","Wq","var_Wq",
                        "L","Lq"])
         
    def solve(self,max_t=10):
        """
        Calculates basic performance results
        
        Returns:
        rho: traffic intensity
        rho_eff: server utilization
        W:  mean time at server
        Wq: mean wait time in queue
        L:  mean number of customers in system
        Lq: mean number of customers in queue
        PK: probability the system is full
        P0: Fraction of time server is idle
        TA: expected number turned away/unit time
        lam_eff: effective arrival rate
        Pn:  matrix of [n,P[N=n]] 
        """
        #traffic intensity
        self.an.traffic_intensity = self.lam/self.mu
        self.an.rho = self.an.traffic_intensity

        #create parameter vector
        n_vec = range(0,self.K+1)
        rho_vec = [self.an.rho] * (self.K+1)
        K_vec = [self.K] * (self.K+1)

        #define Pn function
        Pn = lambda n,rho,K: (1.0 - rho) * rho**n / (1.0 - rho**(K+1.0))

        if self.an.rho != 1.0:
            self.an.Pn = np.column_stack([n_vec,map(Pn,n_vec,rho_vec,K_vec)])
        else:
            self.an.Pn = np.column_stack([n_vec,[1.0/(self.K + 1.0)] * (self.K + 1)])

        #Calculate core results
        self.an.PK = self.an.Pn[self.K,1]
        self.an.lam_eff = self.lam*(1.0-self.an.PK)
        self.an.rho_eff = self.an.lam_eff * self.st
        if self.an.rho != 1.0:
            self.an.Lq = self.an.rho / (1.0 - self.an.rho) - \
                      self.an.rho * (self.K * self.an.rho ** self.K + 1) / \
                      (1.0 - self.an.rho ** (self.K + 1))
            self.an.P0 = (1.0 - self.an.rho) / (1 - self.an.rho**(self.K + 1))
        else:
            self.an.Lq = (self.K*(self.K-1.0))/(2*(self.K+1))
            self.an.P0 = 1.0/(self.K+1)
            
        self.an.L = self.an.Lq + self.an.rho_eff
        self.an.W = self.an.L/self.an.lam_eff
        self.an.Wq = self.an.W - self.st
        self.an.TA = self.lam*self.an.PK

        self.calc_CDF(max_t)

    def calc_Pn(self,max_n):
        pass

    def calc_CDF(self,max_t):
        """
        Calculates basic performance results
        
        Returns:
        CDFWqt: matrix of [t,Wq(t)]
        """
        Q_vec = self.an.Pn[:,1].tolist()
        for n in range(0,self.K):
            Q_vec[n] = Q_vec[n] / (1.0 - self.an.PK)

        #compute Waiting time distribution Wq(t)
        #create parameter vectors
        max_pts = QtsGlobalParameters.max_plot_points
        #create vector of time points (t)
        t_vec = generate_time_vector(max_t,max_pts)
        Wq_tvec = [0.0] * (max_pts+1)

        #compute Wq(t)
        for n in range(0,max_pts+1):
            total = 0.0
            for j in range(1,self.K):
                sumwork = 0.0
                poiss = 1.0
                for i in range(0,j):
                    if i > 0: poiss *= self.mu * t_vec[n] / i
                    sumwork += poiss
                total += Q_vec[j] * sumwork * exp(-self.mu*t_vec[n])
                Wq_tvec[n] = 1.0 - total

        self.CDFWqt = np.column_stack([t_vec,Wq_tvec])

    def simulate(self,max_simulate=1000,max_t=10):
         
        self.simulation_maxNumber = max_simulate

        #set up simulation object
        self.qsim =  BasicSimulator(self.sim,
                          iat_fctn = self.iat_fctn,
                          iat_parms = self.iat_parms,
                          st_fctn = self.st_fctn,
                          st_parms = self.st_parms,
                          K=self.K)

        #run simulation model
        self.qsim.run_simulation(self.simulation_maxNumber,
                                 self.K,max_t)
            
           
if __name__ == "__main__":
    QtsGlobalParameters.simulation_checkpoint=25000
    #QtPyGlobalParameters.TRACING=True
    q1 = M_M_1_K(2.0,0.8,5)
    #q1.st=18.0
    q1.solve()
    q1.print_results()
    q1.calc_CDF(20)
    #q1.plot_Pn(n_incr=2)
    #q1.plot_Pn(show_model_parameters=True, n_incr=2)
    #q1.plot_CDF(show_model_parameters=True)
    q1.simulate(100000,max_t=20)
    q1.print_results("simulation")
    q1.plot_Pn("simulation")
    q1.plot_CDF("simulation")

    x1 = q1.an.Pn[:,0]
    y1 = q1.an.Pn[:,1]
    x2 = q1.sim.Pn[:,0]
    y2 = q1.sim.Pn[:,1]


    side_by_side_barchart(x1,y1,x2,y2,label1="analytic",label2="simulation",n_incr=5)

    
