#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from scipy.interpolate import CubicSpline

#cols=['#FF6E1E','#76777B']
cols=['#FF6E1E','black']

k,p = np.loadtxt('jspectrum.csv',dtype=(float,float), delimiter=",", unpack=True)
csp=CubicSpline(k,p)
kvals = np.logspace(np.log10(k.min()),np.log10(k.max()),100000)
#print kvals
#cs = UnivariateSpline(k,p,k=3,s=0)
#cs = splrep(k,p,k=2,s=0.0)
#csp= PPoly.from_spline(cs,extrapolate=True)

kdown = np.logspace(-4,3,50)
#print kdown
pdown = csp(kdown)

deriv= (csp.derivative(nu=1))(k)
#Setting boundary derivatives to match the finer spline.
#For some reason, the default boundary derivatives caused trouble.
cspdown = CubicSpline(kdown,pdown,bc_type=((1,deriv[0]),(1,deriv[-1])))

print "MAXDIFF: ",np.max(np.abs(np.subtract(csp(kvals),cspdown(kvals))))



#print k
#print len(k)
fig,ax=plt.subplots()
#ax.plot(k,np.multiply(k,p),'-',lw=2)
#ax.plot(k,p,'o',color=cols[1],lw=3, label='400 Samples')
ax.plot(kvals,kvals**4*csp(kvals),'-',color=cols[1],lw=3, label='Cubic Spline')
#ax.plot(kvals,0.9*((np.exp(-20*kvals)+ 0.1*np.exp(-5*kvals))),'-')
#ax.plot(kvals,np.exp(-kvals),'-')
#ax.plot(kdown,pdown,'o',color=cols[0],lw=3,markersize=0, label='100 Samples')
#ax.plot(kvals,cspdown(kvals),'-',color=cols[0],lw=2,markersize=12,label='Cubic Spline (Downsampled)')
fig.suptitle('4x Downsampling',fontsize=14)
ax.set_xlabel(r'$k$')
#ax.set_ylabel(r'$k^2 P(k)$')
ax.set_ylabel(r'$P(k)$')
ax.set_yscale("log")
ax.set_xscale("log")
#ax.set_xlim([1e-4,1.5e-2])
#ax.set_ylim([1e-8,2])

plt.legend(loc='lower left')
plt.show()
