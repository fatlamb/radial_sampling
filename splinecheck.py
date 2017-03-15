#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline, PPoly, splrep,splev


class quadspline:

	def __init__(self,coeff,bkpts):
		self.coeff=coeff
		self.bkpts=bkpts


	def speval(self,x):
		if (x<self.bkpts[0] or x>self.bkpts[-1]):
			print "x out of breakpoint range"

		arg = np.argmin(np.abs(self.bkpts-x))
	
		if (x< self.bkpts[arg]):
			arg -= 1
		
		#x0 = (self.bkpts[arg]+self.bkpts[arg+1])/2
		x0 = self.bkpts[arg] 

		splval = 0
		#splval = self.coeff[2,arg+2]
		for power in range(0,3):
			splval += self.coeff[power,arg]*(x-x0)**(2-power)

		return splval
	
k,p = np.loadtxt('spectrum.csv',dtype=(float,float), delimiter=",", unpack=True)
kvals = np.logspace(np.log10(k.min()),np.log10(k.max()),100000)
print kvals
#cs = UnivariateSpline(k,p,k=3,s=0)
cs = splrep(k,p,k=2,s=0.0)
csp= PPoly.from_spline(cs,extrapolate=True)

print csp.c
print p
print csp.c.shape
print csp.x.shape
print k.shape
qs = quadspline(csp.c,csp.x)
qspeval_vec = np.vectorize(qs.speval)
fig,ax=plt.subplots()

ax.plot(k,p,'o-')
ax.plot(kvals,csp(kvals),'g-')
ax.plot(kvals,qspeval_vec(kvals),'r-')
ax.plot(csp.x,qspeval_vec(csp.x),'ko')
ax.set_yscale("log")
ax.set_xscale("log")
plt.show()
