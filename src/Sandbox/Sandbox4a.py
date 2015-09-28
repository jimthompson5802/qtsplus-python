from cmath import *
import scipy.optimize as opt
import scipy.signal as sig
import numpy as np


lam = 6./5
j = 1
mu = 3./2
k = 2


p = [j*lam]*j + [-k*mu]*k
print "\np=",p
p = np.array(np.poly(p))
#c = (6./5)*3*3
c = ((j*lam)**j)*((k*mu)**k)
print "\np=",p," c=",c

q = np.array(([0]*(len(p)-1))+[c])
print "\nq=",q

r = q - ((-1)**j * p)
print "\nr=",r
#r[len(r)-1] = 0
#print "\nr'=",r

roots = np.roots(r)
print "\nroots=",roots

