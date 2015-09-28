from QtPy.QtPyCommon import *
#Other required imports

class ModelName(BasicQueue):
    """
    <ModelName> Model:  <descrpitor>
    """

    def __init__(self,lam,st):
        """
        Parameters:
        
        lam: mean arrival rate
        st: mean service time
        """
        self.model_name="ModelName"
        BasicQueue.__init__(self)

        self.lam = float(lam)
        self.st = float(st)
        self.c = 1
        self.K = np.Infinity

        
        #set up simulation parameters
        self.iat_fctn = expovariate
        self.st_fctn = expovariate

        #register model specific parameters and results
        self.add_parameter_variable("xxx","Parameter label to(xxx)=%g")

        for obj in [self.an,self.sim]:
            obj.add_result_variable("yyy","results label to print %g")
            


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

    def solve(self):
        """
        Calculates basic performance results

        Returns:
        rho: server utilization
        W:  mean time at server
        Wq: mean wait time in queue
        L:  mean number of customers in system
        Lq: mean number of customers in queue

        """
        pass

    def calc_Pn(self,max_n=10):
        """
        Calculates probability of N=n in the system

        Parameter:
        max_n: maximum system size for plotting

        Returns:
        Pn:  matrix of [n,P[N=n]]
        """
        #define function to calculate P[N=n]
        Pn = lambda n,rho: (1.0 - rho)*rho**n

        #create parameter vectors
        n_vec = range(0,max_n+1)
        #rho_vec = [self.rho]*(max_n+1)

        #create matrxi [n,P[N=n]]
        #self.Pn = np.column_stack([n_vec,map(Pn,n_vec,rho_vec),np.zeros(len(n_vec))])

        #calculate CDF or P[N=n]
        self.Pn[0,2] = self.Pn[0,1]
        for i in xrange(1,len(self.Pn)):
            self.Pn[i,2] = self.Pn[i,1] + self.Pn[i-1,2]

        
        #return self.Pn
 


    def calc_CDF(self,max_t=10):
       """
        Calculates CDF for W(t) and Wq(t)

        Parameter:
        max_t: maximum time for W(t) and Wq(t) plots

        Returns:
        CDFWt: matrix of [t,W(t)]
        CDFWqt: matrix of [t,Wq(t)]
        """
        #create vector of time points (t)
        t_vec = generate_time_vector(max_t,GlobalParameters.max_plot_points)

        #define functions for W(t) and Wq(t)
        #W_t = lambda t,lam,mu: 1.0 - exp(-(mu - lam)*t)
        #Wq_t = lambda t,lam,mu: 1.0 - (lam/mu)*exp(-(mu - lam)*t)

        #create required paramenter vectors
        #lam_vec = [self.lam] * 101
        #mu_vec = [self.mu] * 101

        #compute CDF matrix
        #self.CDF = np.column_stack([t_vec,map(W_t,t_vec,lam_vec,mu_vec),
        #        map(Wq_t,t_vec,lam_vec,mu_vec)])

    def simulate(self,maxNumber=1000,max_n=10,max_t=10):
        
        #set up simulation parameters for arrival and service times
        self.simulation_maxNumber = maxNumber

        #set up simulation object
        self.qsim =  BasicSimulator(self.sim,
                          iat_fctn = self.sim.iat_fctn,
                          iat_parms = self.sim.iat_parms,
                          st_fctn = self.sim.st_fctn,
                          st_parms = self.sim.st_parms,
                          K=self.sim.K)

        #run simulation model
        self.qsim.run_simulation(self.simulation_maxNumber,max_n,max_t)
            
if __name__ == "__main__":
    pass
    #q = ModelName(<>,<>)
    #q.solve()
    #q.print_results()
    #q.calc_Pn(10)
    #q.calc_CDF(10)
    #q.Pn
    #q.CDF
   
    
    
