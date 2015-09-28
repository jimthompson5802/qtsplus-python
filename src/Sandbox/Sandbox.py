from cmath import *
import scipy.optimize as opt
import scipy.signal as sig
import numpy as np



lam = 6./5
mu = 3./2
j= 3
k = 2

roots = []
for n in range(1,k+1):
    f = lambda z: (((j*lam)/(j*lam + k*mu*(1.0-z)))**(float(j)/k)) \
        * exp(complex(0,2.0*n*pi/k))
                                                                                
    if n == k:
        z0 = 0.5
    else:
        z0 = exp(complex(0,2.0*n*pi/float(k)))

    r = opt.fixed_point(f,z0)
    if abs(r.imag) < 1e-8:
        r = r.real
    roots.append(r)            
print "roots:"
for x in roots:
    print x

print "\nConstants:"
C1 = np.array(roots)
C2 = np.array([k*mu]*k) 
C0 = -C1/C2
C = np.prod(C0)
print "C1=",C1
print "C2=",C2
print "C0=",C0
print "C=",C
top1 = np.poly(-C2)
print "top1=",top1
top = C*top1
print "top=",top

bot1 = np.array(roots + [0])
print "bot1=",bot1
bot = np.poly(bot1)
print "bot=",bot

print "partial fraction expansion\n"
for x in sig.residue(top,bot):
    print x



