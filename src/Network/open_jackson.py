import numpy as np
import math as mt

def open_jackson(gamma_vec, mu_vec, c_vec, R):
    """
    Open Jackson Network:  N M/M/c network routing matrix
    
    Parameters:
        gamma_vec: vector of external arrival rates
        mu_vec: vector of service rates
        c_vec: vector of number of servers
        R: routing matrix

    Results:
        lam_vec: node arrival rate vector
        rho_vec: server utilization vector
        W:  soujourn time through system
        L:  mean number of customers in system
        W_vec: vector of node wait times
        Wq_vec: vector of node queue wait times
        L_vec: vector of mean number of customers at node
        Lq_vec: vector of mean number of customers in node queue
        
    """

    #solve traffic equations
    rows,cols=R.shape
    lam_vec = np.linalg.solve(-(R.T-np.identity(rows)),gamma_vec)

    #calculate core results
    rho_vec = lam_vec/(mu_vec*c_vec)

    #build result dictionary
    ans = {"lam_vec":lam_vec,"rho_vec":rho_vec}

    print P_Wait(lam_vec,mu_vec,c_vec)

    return ans


def P_Wait(lam, mu, c):
    
    P0 = 1.0
    X = 1.0

    for n in range(1,c):
        X = X*(lam/mu)/n
        P0=P0+X
        

    X = X * (lam / mu) / c
    X = X / (1 - (lam / (mu * c)))
    
    P0 = 1 / (P0 + X)
    
    return X * P0





if __name__ == "__main__":
    gamma_vec=np.array([0.,0.,0.,4.])
    mu_vec=np.array([25.,100./3,50./3,20.])
    c_vec=np.array([2,1,1,1])
    R=np.array([[0.,0.5,0.5,0.],
        [1.,0.,0.,0.],
        [0.6,0.,0.,0.],
        [1.0,0.,0.,0.]])
    print open_jackson(gamma_vec,mu_vec,c_vec,R)
