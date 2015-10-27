import numpy as np
import matplotlib.pyplot as plt
import SimPy.Simulation as sp
from random import expovariate

###
# Class to hold global parameters that affect run-time
###
class QtsGlobalParameters:
    """
    Container for the following parameters
    epsilon:  tolerance to determine convergence for numerical routines
    max_iteration:  number of iterations for iterative numerical routines
    max_plot_points: number of plot points for time distribution plots
    TRACING: flag to turn on/off tracing during simulation
    simulation_checkpoint: frequency to print simulation progress message
    simulation_sampling_interval: simulated time frequency for sampling
        system size and queueu size metrics
    """
    #convergence tolerance for iterative methods
    epsilon = 1e-8
    #maximum number of interations allowed for iterative methods
    max_iteration = 100

    #number of plot points for CDF plotting
    max_plot_points = 100

    #simulation tracing flag
    TRACING = False               

    #Print checkpoint message every completion, 0 to turn off checkpoint
    simulation_checkpoint = 1000  

    #mean time for an exponential time interval to sample system size
    simulation_sampling_interval = 10

    
###
#instantiate Global Parameter object
###
QtPyGlobalParameters = QtsGlobalParameters()


###
# Container class for model results
###
class ModelResults():

    def __init__(self):
        # dictionary containing labels for each result variable
        self.result_labels = {}

    def add_result_variable(self,name,print_label):
        self.result_labels[name] = print_label
        self.__dict__[name] = None

    def print_results(self, result_sequence):
        for var in result_sequence:
            if self.__dict__.get(var) != None:
                print ("    " + self.result_labels[var])% self.__dict__.get(var)


###
# Base class for all queueing models
###
class BasicQueue():
    """
    Parent Class for non-networked queueing models
    """
    model_name = ""
    def __init__(self):
        self.gp = QtsGlobalParameters()
        # dictionary containing labels for each parameter
        self.parameter_labels = {}

        #create parameter labels     
        self.add_parameter_variable("lam","Arrival rate(lam)=%g")
        self.add_parameter_variable("iat","Mean inter-arrival time=%g")
        self.add_parameter_variable("st","Mean service time(st)=%g")
        self.add_parameter_variable("mu","Service rate(mu)=%g")
        self.add_parameter_variable("iat_fctn","Inter-arrival time function: %s")                                       
        self.add_parameter_variable("iat_parms","Inter-arrival time function parameters: %s")                                       
        self.add_parameter_variable("st_fctn","Service time function: %s")                                       
        self.add_parameter_variable("st_parms","Service time function parameters: %s")                                       
        self.add_parameter_variable("c","number of servers(c)=%g")
        self.add_parameter_variable("K","Maximum capacity of system(K)=%d")

        #establish area for analytic and simulation model results and related labels
        self.an = ModelResults()    #container for analytic model results
        self.sim = ModelResults()   #container for simulation model results
        for obj in [self.an, self.sim]:                                       
            obj.add_result_variable("end_simulation_time","Simulation ended at %g")
            obj.add_result_variable("transactions_arrived","Number of transaction arrivals(simulation): %g")
            obj.add_result_variable("transactions_completed","Number of transaction completed(simulation): %g")
            obj.add_result_variable("rho","Server utilization(rho)=%g")
            obj.add_result_variable("rho_eff","Effective server utilization(rho_eff)=%g")
            obj.add_result_variable("avg_iat","Average of generated inter-arrival time(avg_iat)=%g")
            obj.add_result_variable("var_iat","Variance of generated inter-arrival time(var_iat)=%g")
            obj.add_result_variable("avg_st","Average generated service time(avg_st)=%g")
            obj.add_result_variable("var_st","Variance of generated service time(var_st)=%g")
            obj.add_result_variable("W","Mean waiting time(W)=%g")
            obj.add_result_variable("var_W","Variance waiting time(varW)=%g")
            obj.add_result_variable("Wq","Mean waiting time in queue(Wq)=%g")
            obj.add_result_variable("var_Wq","Variance waiting time in queue(varWq)=%g")
            obj.add_result_variable("L","Mean number of customers in the system(L)=%g")
            obj.add_result_variable("Lq","Mean number of customers in the queue(Lq)=%g")
            obj.add_result_variable("P0","Fraction of time server is idle(P0)=%g")
            obj.add_result_variable("PK","Probability that the systems is full(PK)=%g")

        
    def solve(self):
        pass

    def __print_value(self,lab,v):
        v1 = self.__dict__.get(v)
        if v1 != None:
            print ("    " + lab) % v1

    def __generate_model_parmater_string(self):
        parm_list = ["lam","st","c","K"]
        parm_string = ""
        for p in parm_list:
            v = self.__dict__.get(p)
            if v != None:
                if v != np.Infinity:
                    if len(parm_string) > 0:
                        parm_string += (", " + p + "=" + ("%g" % v))
                    else:
                        parm_string = p + "=" + ("%g" % v)
        return parm_string

    def add_parameter_variable(self,name,print_label):
        self.parameter_labels[name] = print_label

    def print_parameters(self,parameter_sequence):
        for var in parameter_sequence:
            v = self.__dict__.get(var)
            if v != None and v != np.Infinity:
                print ("    " + self.parameter_labels[var])% v

    def print_results(self,type="analytic"):
        """
        Print model parameters and available analytic results
        """
        if type == "analytic":
            print "\nModel:  %s" % self.model_name
            print "Parameters:"
            self.__print_model_specific_parameters__(type)
            print "Results:"
            self.__print_model_specific_results__(type)
        
        elif type == "simulation":
            print "\nSimulation of Model:  %s" % self.model_name
            print "Parameters:"
            self.__print_model_specific_parameters__(type)
            print "Results:"
            self.__print_model_specific_results__(type)

        else:
            raise ParameterError("Invalid value for 'type' parameter: %s" % type)



    def plot_Pn(self,type="analytic",max_n=None,n_incr=1,bar_width = 0.35,show_model_parameters=False):
        if type == "analytic":
            obj = self.an
        elif type == "simulation":
            obj = self.sim
        else:
            raise ParameterError("Invalid value for 'type' parameter: %s" % type)
        
        if max_n != None:
            self.calc_Pn(max_n)
        else:
            max_n = max(obj.Pn[:,0])

        x = obj.Pn[0:max_n+1,0]
        y = obj.Pn[0:max_n+1,1]
        
        width = bar_width
        p1 = plt.bar(x, y, width, color='r')

        #setup labels for x-coordinates
        nlab = np.arange(min(x),max(x)+1,n_incr)

        #set up title
        if show_model_parameters:
            subtitle = "\n" + self.__generate_model_parmater_string()
        else:
            subtitle = ""
        
        plt.ylabel('P[N=n]')
        plt.xlabel('N')
        plt.title('System Size Distribution for ' + self.model_name + subtitle)
        plt.xticks(nlab+width/2., map(int,nlab) )

        plt.show()

    def plot_CDF(self,type="analytic",show_model_parameters=False,legend_loc="best"):
        if type == "analytic":
            obj = self.an
            col = 1
        elif type == "simulation":
            obj = self.sim
            col = 2
        else:
            raise ParameterError("Invalid value for 'type' parameter: %s" % type)
        
        CDFWt = obj.__dict__.get("CDFWt")
        if CDFWt != None:
            x1 = CDFWt[:,0]
            y1 = CDFWt[:,col]
            plt.plot(x1,y1,"g-",label="W(t)",linewidth=1.5)

        CDFWqt = obj.__dict__.get("CDFWqt")
        if CDFWqt != None:
            x2 = CDFWqt[:,0]
            y2 = CDFWqt[:,col]
            plt.plot(x2,y2,"r--",label="Wq(t)",linewidth=1.5)

        #set up title
        if show_model_parameters:
            subtitle = "\n" + self.__generate_model_parmater_string()
        else:
            subtitle = ""

        plt.title('CDF Plot for ' + self.model_name + subtitle)
        plt.legend(loc=legend_loc)
        plt.show()


    def simulate(self,maxNumber=1000):
        pass

    def reset_simulation(self):
        del(self.qsim)
        
        

###
# Wrapper for SimPy
###
class BasicSimulator(sp.Simulation):
    """
    Simulation class
    Parameters:
    iat_fctn: function object for generating inter-arrival times
    iat_params:  list of parameters used by iat_fctn
    st_fctn:  function object for generating service times
    st_parms: list of parameters used by st_fctn
    c:  number of servers
    K: capacity of queueing system (in service and in queue)
    """   
    def __init__(self,sim,iat_fctn,iat_parms,st_fctn,st_parms,c=1,K=np.Infinity):
        self.iat_fctn = iat_fctn
        self.iat_parms = iat_parms
        self.st_fctn = st_fctn
        self.st_parms = st_parms
        self.c = c
        self.K = K
        self.sim = sim       #simulation model results container

    def run_simulation(self,max_simulate=1000,max_n=10,max_t=10):
        """
        Starts the simulation run
        max_simulate: maximun number of transactions to simulate
        max_n: number of bins for the Pn matrix
        max_t: upper time limit for collecting W(t) and Wq(t) statistics
        """
        self.maxNumber = max_simulate
        self.initialize()

        print "\nstarting simulation for " + `max_simulate` + " transactions"

        #Server Definition        
        self.Server=sp.Resource(capacity=self.c,name='Server',
                             monitored=False,
                             sim=self)

        self.__setup_monitors(max_n,max_t)

        #run the simulation
        TA = TransactionArrival("Trans",sim=self)
        self.activate(TA,TA.arrival(iat_fctn=self.iat_fctn,
                                    iat_parms=self.iat_parms,
                                    st_fctn=self.st_fctn,
                                    st_parms=self.st_parms,
                                    K=self.K))
        
        MS = MonitorSystemSize("MonSize",sim=self)
        self.activate(MS,MS.record_system_size())
        
        self.simulate(until=sp.infinity)

        self.__record_performance_metrics()

        print "simulation completed"

    def __setup_monitors(self,max_n,max_t):
        
        #Monitors for key metrics                        
        self.mon_Pn=sp.Tally(sim=self)
        self.mon_Pn.setHistogram("Pn",0,max_n+1,max_n+1)
        self.mon_Qn=sp.Tally(sim=self)
        self.mon_Qn.setHistogram("Qn",0,max_n+1,max_n+1)
        #self.mon_L = sp.Tally(sim=self)       #number in system
        #self.mon_Lq = sp.Tally(sim=self)      #number in queue
        self.mon_W = sp.Tally(sim=self)       #time in system
        nbins = QtPyGlobalParameters.max_plot_points
        self.mon_W.setHistogram("W(t)",0,max_t,nbins)
        self.mon_Wq = sp.Tally(sim=self)      #time in queue
        self.mon_Wq.setHistogram("Wq(t)",0,max_t,nbins)
        self.mon_Server_Busy = sp.Tally(sim=self)  #time the server is busy
        self.mon_iat = sp.Tally(sim=self)     #inter-arrival time
        self.mon_arrivals = sp.Tally(sim=self)  #count of arrivals
        self.mon_balks = sp.Tally(sim=self)   #count of transactions that balk
            

    def __record_performance_metrics(self):
        
        #record key performance metrics
        self.sim.end_simulation_time = self.now()
        self.sim.transactions_completed = self.mon_W.count()
        self.sim.transactions_arrived = self.mon_arrivals.total()
        self.sim.rho = self.mon_Server_Busy.total()/(self.c * self.now())
        self.sim.W = self.mon_W.mean()
        self.sim.var_W = self.mon_W.var()
        self.sim.Wq = self.mon_Wq.mean()
        self.sim.var_Wq = self.mon_Wq.var()
        #self.sim.L = self.mon_L.timeAverage()
        self.sim.L = self.mon_Pn.mean()
        #self.sim.Lq = self.mon_Lq.timeAverage()
        self.sim.Lq = self.mon_Qn.mean()
        self.sim.avg_iat = self.mon_iat.mean()
        self.sim.var_iat = self.mon_iat.var()
        self.sim.avg_st = self.mon_Server_Busy.mean()
        self.sim.var_st = self.mon_Server_Busy.var()

        #calculate Pn matrix
        len_histogram = len(self.mon_Pn.getHistogram())
        self.Pn = np.array(self.mon_Pn.getHistogram()[1:len_histogram-1])
        self.Pn[:,1] = self.Pn[:,1]/float(self.mon_Pn.count())
        self.Pn = np.column_stack([self.Pn,np.zeros(len(self.Pn))])
        self.Pn[0,2] = self.Pn[0,1]
        for i in xrange(1,len(self.Pn)):
            self.Pn[i,2] = self.Pn[i,1] + self.Pn[i-1,2]
        self.sim.Pn = self.Pn  #copy to simulation model results area
        del(self.Pn)

        #calculate CDFWt
        len_histogram = len(self.mon_W.getHistogram())
        self.CDFWt = np.array(self.mon_W.getHistogram()[1:len_histogram-1])
        self.CDFWt[:,1] = self.CDFWt[:,1] / float(self.mon_W.count())
        self.CDFWt = np.column_stack([self.CDFWt,np.zeros(len(self.CDFWt))])
        self.CDFWt[0,2] = self.CDFWt[0,1]
        for i in xrange(1,len(self.CDFWt)):
            self.CDFWt[i,2] = self.CDFWt[i,1] + self.CDFWt[i-1,2]
        self.sim.CDFWt = self.CDFWt     #copy to simulation model results area
        del(self.CDFWt)
        
        #calculate CDFWqt
        len_histogram = len(self.mon_Wq.getHistogram())
        self.CDFWqt = np.array(self.mon_Wq.getHistogram()[1:len_histogram-1])
        self.CDFWqt[:,1] = self.CDFWqt[:,1] / float(self.mon_Wq.count())
        self.CDFWqt = np.column_stack([self.CDFWqt,np.zeros(len(self.CDFWqt))])
        self.CDFWqt[0,2] = self.CDFWqt[0,1]
        for i in xrange(1,len(self.CDFWqt)):
            self.CDFWqt[i,2] = self.CDFWqt[i,1] + self.CDFWqt[i-1,2]
        self.sim.CDFWqt = self.CDFWqt
        del(self.CDFWqt)
        


class TransactionArrival(sp.Process):
    """
    Defines arrival process for transaction
    """
    def arrival(self,iat_fctn,iat_parms,st_fctn,st_parms,K):
        i = 0
        while True:
            #create a transactions
            i += 1
            T = TransactionWork("Transid " + `i`,sim=self.sim)

            #Insert newly created transaction into server
            self.sim.activate(T,T.run(st_fctn,st_parms,K),delay=0)

            #wait for next transaction arrival
            iat = map(iat_fctn,*iat_parms)[0]
            self.sim.mon_iat.observe(iat)
            yield sp.hold,self,iat

class MonitorSystemSize(sp.Process):
    """
    Process to monitor number of transactions in the system during simulation
    This monitor is based on Poisson Arrivals See Time Averages(PASTA) principle
    """
    def record_system_size(self):
        t = float(QtPyGlobalParameters.simulation_sampling_interval)
        while True:
            #collect sample of number of transactions in the system
            self.sim.mon_Pn.observe(len(self.sim.Server.waitQ)
                                + len(self.sim.Server.activeQ))
            self.sim.mon_Qn.observe(len(self.sim.Server.waitQ))

            #wait for next sampling time
            iat = expovariate(1/t)
            yield sp.hold,self,iat
            
class TransactionWork(sp.Process):
    """
    defines what the a transaction performs at Server
    """
    def run(self,st_fctn,st_parms,K):
        
        self.__trace("arrived   ")
        
        #record arrival time of new transaction
        arrival_time = self.sim.now()
        self.sim.mon_arrivals.observe(1)

        #deterimine if transaction can enter system or not
        if (K == np.Infinity or
            (len(self.sim.Server.waitQ)+len(self.sim.Server.activeQ))) < K:
            #transaction can enter system
            #record number in system, including the new arrived transaction ("+1")
            #self.sim.mon_L.observe(len(self.sim.Server.waitQ)+len(self.sim.Server.activeQ)+1)
            
            #record number in queue
            #self.sim.mon_Lq.observe(len(self.sim.Server.waitQ))   

            #request server, if not available enter server queue
            yield sp.request,self,self.sim.Server

            #record time in queue
            self.sim.mon_Wq.observe(self.sim.now()-arrival_time)  

            #calculate time that server will be busy and record the server busy time
            service_time = map(st_fctn,*st_parms)[0]  
            self.sim.mon_Server_Busy.observe(service_time)  

            #transaction using server
            self.__trace("using server   ")
            yield sp.hold,self,service_time        
            
            #transaction completed processing at server
            yield sp.release,self,self.sim.Server          

            #capture statistics at completion of transaction
            #self.sim.mon_L.observe(len(self.sim.Server.waitQ)+len(self.sim.Server.activeQ))
            self.sim.mon_W.observe(self.sim.now()-arrival_time)  

            #deterimine if checkpoint message should be written
            if ((QtPyGlobalParameters.simulation_checkpoint != 0) and
              (self.sim.mon_W.count() % QtPyGlobalParameters.simulation_checkpoint) == 0):
                print "Completed %d transactions" % self.sim.mon_W.count()

            #at end of simulation?
            self.__trace("leaving server   ")
            if self.sim.mon_W.count() >= self.sim.maxNumber:
                self.sim.stopSimulation()


        else:
            #transaction unable to enter system, record that we had a balk
            self.__trace("balk   ")
            self.sim.mon_balks.observe(1)


    def __trace(self,message):
        FMT="%7.4f %6s %10s (%2d)"
        if QtPyGlobalParameters.TRACING:
            print FMT%(self.sim.now(),self.name,message,self.sim.mon_L._last_observation)

###
# Common functions
###
def side_by_side_barchart(x1,y1,x2,y2,
                          color1="r",color2="g",
                          label1="",label2="",
                          max_n=None,n_incr=1,bar_width = 0.35,title=""):

    if max_n != None:
        max_n = max(x1)
    
    p1 = plt.bar(x1, y1, width=bar_width, color=color1,label=label1)
    p2 = plt.bar(x2+bar_width, y2, width=bar_width, color=color2,label=label2)

    #setup labels for x-coordinates
    nlab = np.arange(min(x1),max(x1)+1,n_incr)
   
    plt.ylabel('P[N=n]')
    plt.xlabel('N')
    plt.legend(loc="best")
    plt.title(title)
    plt.xticks(nlab+bar_width, map(int,nlab) )
    plt.show()

def side_by_side_plot(x1,y1,x2,y2,
                      color1="r",color2="g",
                      label1="", label2="",
                      max_t=None,t_incr=1,
                      title=""):

    if max_t != None:
        max_t = max(x1)
        
    p1 = plt.plot(x1,y1,color=color1,label=label1)
    p2 = plt.plot(x2,y2,color=color2,label=label2)

    #set up labels for x-coordinates
    nlab = np.arange(min(x1),max(x1),t_incr)

    plt.xlabel("t")
    plt.legend(loc="best")
    plt.title(title)
    plt.xticks(nlab, map(int,nlab))
    plt.show()
    


def erlang_k(mean,k):
    """
    Puedo-random number for Erlang-k distribution
    Parameters:
    mean:  mean for the Erlang-k distribution
    k: Shape parameter
    """
    ans = 0.0
    for i in range(1,k+1):
        mt = mean/float(k)
        ans += expovariate(1.0/mt)
    return ans

def generate_time_vector(max_t, max_pts):
    incr = max_t/float(max_pts)
    n_vec = range(0,max_pts+1)
    t_vec = [incr] * (max_pts+1)
    t_vec = map(lambda n,t: n*t,n_vec,t_vec)
    return t_vec
    
def ErlangB(c,r):
    """
    Computes Erlang-B formula
    
    """
    if c == 0:
        ans = 1.0
    else:
        ErlBTemp = ErlangB(c-1,r)
        ans = (r * ErlBTemp) / (c + r * ErlBTemp)
    return ans

def ErlangC(c,r):
    """
    Compute Erlang-C formula
    """
    ErlBTemp = ErlangB(c,r)
    return (c * ErlBTemp) / (c - r + r * ErlBTemp)

###
# User defined Exceptions
###
class ServerUnstable(Exception):
    def __init__(self,arg = None):
        self.arg = arg
    def __str__(self):
        return str(self.arg)

class ParameterError(Exception):
    def __init__(self,arg=None):
        self.arg = arg
    def __str__(self):
        return str(self.arg)
    
