import QtPy
QtPy.QtPyGlobalParameters.simulation_checkpoint = 0
MaxSimulate = 100000
for j in range(1,11):
    for k in range(1,11):
        q = QtPy.Ej_Ek_1(6./5,j,2./3,k)
        print "j=",j,", k=",k
        q.solve()
        #q.simulate(MaxSimulate)
        #print "xxx",j,k,q.an.rho,q.sim.rho,q.an.W,q.sim.W,q.an.Wq,q.sim.Wq,q.an.L,q.sim.L,q.an.Lq,q.sim.Lq
        
        

