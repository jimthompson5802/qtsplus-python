from cmath import *
import scipy.optimize as opt
import scipy.signal as sig
import numpy as np



lam = 6./5
j= 1
mu = 3./2
k= 2

p = []
for n in range(j+k+1):
    if n == j+k:
        p.append(k*mu)
    elif n == k:
        p.append(-(j*lam+k*mu))
    elif n == 0:
        p.append(j*lam)
    else:
        p.append(0)
print "p=",p
print p.reverse()
print "p=",p

roots = np.roots(p)
print "roots=",roots
        
z = []
for x in roots:
    if x.real < 0:
        z.append(x)
print "z=",z

bot = np.poly(z + [0])
print "bot=",bot

C1 = np.array([k*mu]*k)
C2 = -np.array(z)/C1
C = np.prod(C2)
print C
print "top1=",np.poly([-k*mu]*k)

top = C* np.poly([-k*mu]*k)
print "top=",top


pfe = sig.residue(top,bot)
r = pfe[0]
p = pfe[1]


print "residue=",r
print "poles=",p

Wqt = lambda t: np.sum(r*np.exp(p*t))

print "Wqt=",Wqt(4.0)



