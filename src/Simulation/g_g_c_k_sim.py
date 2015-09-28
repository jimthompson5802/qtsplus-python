from QtPy.QtPyCommon import *
#Other required imports
import numpy as np
import random

class G_G_c_K_Sim(BasicQueue):
    """
    G/G/c Simulation Model:  General inter-arrival and service times, multiple servers
    Parameters:
    iat_fctn:  function object to generate inter-arrival times
    iat_parms:  list containing parameters for the iat_fctn
    st_fctn:  function object to generate the service times
    st_parms:  list containing parameters for the st_fctn
    c:  number of servers (Default = 1)
    K:  maximum number of transactions allowed in system (default=infinity)
    """

    def __init__(self,iat_fctn,iat_parms,st_fctn,st_parms,c=1,K=np.Infinity):
        """
        Parameters:
        
        lam: mean arrival rate
        st: mean service time
        """
        self.model_name="G/G/c/K"
        self.iat_fctn = iat_fctn
        self.iat_parms = iat_parms
        self.st_fctn = st_fctn
        self.st_parms = st_parms
        self.c = c
        self.K = K


    def simulate(self,max_simulate=1000,max_n=10,max_t=10):
        
        #set up simulation parameters for arrival and service times
        self.simulation_maxNumber = max_simulate

        #set up simulation object
        self.qsim =  BasicSimulator(iat_fctn = self.iat_fctn,
                          iat_parms = self.iat_parms,
                          st_fctn = self.st_fctn,
                          st_parms = self.st_parms,
                          c=self.c,
                          K=self.K)

        #run simulation model
        self.qsim.run_simulation(self.simulation_maxNumber,
                                 max_n,max_t)

if __name__ == "__main__":
    QtPyGlobalParameters.TRACING=False
    QtPyGlobalParameters.simulation_checkpoint=25000
    D = lambda x: x
    q = G_G_c_K_Sim(
        expovariate,
        [[1/3.0]],
        expovariate,
        [[1/0.2]],
        c=1,
        K=10)
    q.simulate(100000,10,10)
    q.print_simulation_results()
    q.plot_simulation_Pn()
    q.plot_simulation_CDF()
    
    
