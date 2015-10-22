import QtPyCommon as G
#Other required imports
import cmath as cm
import numpy as np
import scipy.optimize as opt
import scipy.signal as sig


class Ej_Ek_1(G.BasicQueue):
    """
    Ej/Ek/1 Model:  Erlang input to single-server queue, with Erlang service,
        unlimited capacity and FIFO discipline.
        See GSTH, Section 6.2.1 for details of analytic solution.
    """

    def __init__(self,lam,j,st,k):
        """
        Parameters:
        
        lam: mean arrival rate
        j: Erlang shape paramter for arrival process
        st: mean service time
        k: Erlang shape parameter for service process
        """
        self.model_name="Ej/Ek/1"
        G.BasicQueue.__init__(self)

        #set up simulation parameters
        self.iat_fctn = G.erlang_k
        self.st_fctn = G.erlang_k
        self.iat_parms = [[],[]]
        self.st_parms = [[],[]]

        #capture paramters
        self.lam = float(lam)
        self.j = j
        self.st = float(st)
        self.k = k
        self.c = 1
        self.K = np.Infinity

        #register model specific parameter labels
        self.add_parameter_variable("j","Shape parameter for arrival distribution(j)=%d")
        self.add_parameter_variable("k","Shape parameter for service distribution(k)=%d")

        
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
            self.__dict__["iat_parms"][0] = [1/float(value)]
        elif item == "j":
            self.__dict__["iat_parms"][1] = [value]
        elif item == "st":
            self.__dict__["mu"] = 1/float(value)
            self.__dict__["st_parms"][0] = [float(value)]
        elif item == "k":
            self.__dict__["st_parms"][1] = [value]

    def __print_model_specific_parameters__(self,type):
        #print model specific parameters
        if type == "analytic":
            self.print_parameters(["lam","j","st","k"])
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
        roots: roots of the characteristic equation

        """
        #initialize empty list to hold roots
        self.an.roots =[]
        lam = self.lam
        mu = self.mu
        k = self.k
        j = self.j
        
        #check for stable queueing system
        self.an.rho = lam/mu
        if self.an.rho >= 1.0:
            raise G.ServerUnstable("rho greater or equal to 1.0, rho=%g" % self.an.rho)

        ###
        #calculate roots of characteristic equation
        #Based on discussion of GEj/GEk/1 queue in GSTH, Section 6.2.1
        ###
        p = np.array(np.poly(([j*lam]*j)+([-k*mu]*k)))
        c = ((j*lam)**j)*((k*mu)**k)
        q = np.array(([0]*(len(p)-1))+[c])

        #(-1**j) required to compensate for np.poly() behaviour of normalizing
        #leading coefficient to 1
        r = q - (((-1)**j) * p)
        print "r[constant]=",r[len(r)-1]
        r[len(r)-1] = 0.          #debugging code

        #compute roots of the charateristic equation
        roots = np.roots(r)

        #extract roots with negative real parts
        z = []
        for x in roots:
            if x.real < 0:
                z.append(x)

        #create bot half of rational polynomial
        bot = np.array(np.poly(z + [0]))

        #create numerator of polynomial
        C1 = np.array([k*mu]*k)
        C2 = -np.array(z)/C1
        C = np.prod(C2)
        top = C * np.poly([-k*mu]*k)

        #create partial fraction expansion of characteristic equation
        pfe = sig.residue(top,bot)
        self.an.resid = pfe[0]         #residue of partial fraction expansion
        self.an.poles = pfe[1]         #poles from partial fraction expansion

        #compuate Wq
        self.an.Wq = 0.
        for n in range(len(self.an.resid)):
            if self.an.poles[n] != 0.0:      #is this "not equal" to zero
                self.an.Wq += self.an.resid[n]/self.an.poles[n]
        self.an.Wq = self.an.Wq.real
                    
        #compute other performance measures
        self.an.W = self.an.Wq + self.st

        #use Little's Law to compute L and Lq
        self.an.Lq = self.lam * self.an.Wq
        self.an.L = self.lam * self.an.W

    #def calc_Pn(self,max_n=10):
        """
       Calculates probability of N=n in the system

        Parameter:
        max_n: maximum system size for plotting

        Returns:
        Pn:  matrix of [n,P[N=n]]

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
        """ 


    def calc_CDF(self,max_t=10):
        """
        Calculates CDF for W(t) and Wq(t)

        Parameter:
        max_t: maximum time for W(t) and Wq(t) plots

        Returns:
        CDFWqt: matrix of [t,Wq(t)]
        """
        #create vector of time points (t)
        max_pts = G.GlobalParameters.max_plot_points
        t_vec = G.generate_time_vector(max_t,max_pts)
        self.an.CDFWqt = np.zeros(max_pts+1)
        r = self.an.resid
        p = self.an.poles
        for n in range(max_pts+1):
            self.an.CDFWqt[n] = (np.sum(r*np.exp(p*t_vec[n]))).real

        self.an.CDFWqt = np.column_stack([t_vec,self.an.CDFWqt])
            
    def simulate(self,maxNumber=1000,max_n=10,max_t=10):
        
        #set up simulation parameters for arrival and service times
        self.simulation_maxNumber = maxNumber

        #set up simulation object
        self.qsim =  G.BasicSimulator(self.sim,
                          iat_fctn = self.iat_fctn,
                          iat_parms = self.iat_parms,
                          st_fctn = self.st_fctn,
                          st_parms = self.st_parms)

        #run simulation model
        self.qsim.run_simulation(self.simulation_maxNumber,max_n,max_t)
            
if __name__ == "__main__":
    G.QtPyGlobalParameters.simulation_checkpoint = 50000
    q = Ej_Ek_1(6./5,5,2./3,3)   #example 6.3
    #q = Ej_Ek_1(2,3,0.3,2)      #Example 6.1
    q.solve()
    q.print_results()
    #q.calc_Pn(10)
    q.calc_CDF(10)
    q.plot_CDF()
    #q.Pn
    #q.CDF
    q.simulate(200000,max_t=5.0)
    q.print_results("simulation")
    G.side_by_side_plot(q.an.CDFWqt[:,0],q.an.CDFWqt[:,1],
                      q.sim.CDFWqt[:,0],q.sim.CDFWqt[:,2],
                      label1="analytic",label2="simulation")
    
    
