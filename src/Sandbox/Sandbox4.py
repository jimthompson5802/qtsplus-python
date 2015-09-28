from cmath import *
import scipy.optimize as opt
import scipy.signal as sig
import numpy as np


lam = 6./5
j = 1
mu = 3./2
k = 21


p = np.array(np.poly(([j*lam]*j)+([-k*mu]*k)))
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

print "\nWqt(2.5)=",Wqt(2.5)
print "\nWqt(3.0)=",Wqt(3.0)
print "\nWqt(4.0)=",Wqt(4.0)

Wq = 0.
for n in range(len(r)):
    if p[n] != 0:
        Wq += r[n]/p[n]

print "\nWq=",Wq
