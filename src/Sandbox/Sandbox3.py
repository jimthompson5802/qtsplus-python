from cmath import *
import scipy.optimize as opt
import scipy.signal as sig
import numpy as np


lam = 6./5
j= 1
mu = 3./2
k= 2

p1 = [lam*j]*j
p1 = [0] + p1
p2 = [-mu*k]*k
p3 = np.array([0]*(j+k+1))
p3[j+k] = (j*lam)**j *  (k*mu)**k
print "\np1+p2=",p1+p2
print "\np3=",p3

bot = np.array(np.poly(p1 + p2))
print "\n1st=",bot

bot = p3 - bot
print "\n2nd=",bot

roots = np.roots(bot)
print "\nroots=",roots

z = []
for x in roots:
    if x.real < 0:
        z.append(x)
print "\nz=",z

bot = np.array(np.poly(z + [0]))
print "\nbot=",bot

C1 = np.array([k*mu]*k)
print "\n-z=",-np.array(z)
C2 = -np.array(z)/C1
print "\nC1=",C1," C2=",C2
C = np.prod(C2)
print "\nC=",C
print "\ntop1=",np.poly([-k*mu]*k)

top = C* np.poly([-k*mu]*k)
print "\ntop=",top


pfe = sig.residue(top,bot)
r = pfe[0]
p = pfe[1]


print "\nresidue=",r
print "\npoles=",p

Wqt = lambda t: np.sum(r*np.exp(p*t))

print "\nWqt(4.0)=",Wqt(4.0)

