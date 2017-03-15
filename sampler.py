#!/usr/bin/python

from scipy.interpolate import PPoly, splrep
from Crawlers import KCrawlerDict
from Crawlers import HCrawlerDict
from Crawlers import XYCrawlerDict
import numpy as np

k,p = np.loadtxt('spectrum.csv',dtype=(float,float), delimiter=",", unpack=True)

#Fit a quadratic spline to the data.
qs = splrep(k,p,k=2,s=0.0)

#Convert this spline to a piecewise quadratic
qsp= PPoly.from_spline(qs,extrapolate=True)


#Generate a dictionary of XYn. This will fill as Xn/Yn(x) are requested.
#These calls are also idempotent, so no values are calculated twice.  
xyca = XYCrawlerDict()

#Generate a dictionary of Kln. 
hca = HCrawlerDict(xyca)


def calc_var(l,r_alpha,local_qsp,local_ca):

	#Extract breakpoints from the piecewise polynomial.
	bkpts = local_qsp.x
	
	#Extract polynomial coefficients.
	coeffs = local_qsp.c
	
	npoints = len(bkpts)
	delta_h = np.zeros((3,npoints-1))
	hlast = np.zeros(3)
	hcurrent = np.zeros(3)
	
	for i in range(npoints):
		for n in range(3):
			hcurrent = local_ca.GetEntry((bkpts[i],r_alpha),(l,n))	
			if i>0:
				delta_h[n,i-1] = hcurrent-hlast
			hlast = hcurrent
		
	
	accum = 0
	for i in range(npoints-1):
		accum += (coeffs[0,i] - coeffs[1,i]*bkpts[i] + coeffs[2,i]*bkpts[i]**2)*delta_h[0,i]
		accum += (coeffs[1,i] - 2*coeffs[2,i]*bkpts[i]**2)*delta_h[1,i]
		accum += coeffs[2,i]*delta_h[2,i]
	
	var = 16.0*np.pi**2/(r_alpha**l * r_beta**l)*accum
	return var



#Generate a dictionary of Kln. 
kca = KCrawlerDict(xyca)


def calc_covar(l,r_alpha,r_beta,local_qsp,local_ca):

	#Extract breakpoints from the piecewise polynomial.
	bkpts = local_qsp.x
	
	#Extract polynomial coefficients.
	coeffs = local_qsp.c
	
	npoints = len(bkpts)
	delta_k = np.zeros((3,npoints-1))
	klast = np.zeros(3)
	kcurrent = np.zeros(3)
	
	for i in range(npoints):
		for n in range(3):
			kcurrent = local_ca.GetEntry((bkpts[i],r_alpha,r_beta),(l,n))	
			if i>0:
				delta_k[n,i-1] = kcurrent-klast
			klast = kcurrent
		
	
	accum = 0
	for i in range(npoints-1):
		accum += (coeffs[0,i] - coeffs[1,i]*bkpts[i] + coeffs[2,i]*bkpts[i]**2)*delta_k[0,i]
		accum += (coeffs[1,i] - 2*coeffs[2,i]*bkpts[i]**2)*delta_k[1,i]
		accum += coeffs[2,i]*delta_k[2,i]
	
	covar = 16.0*np.pi**2/(r_alpha**l * r_beta**l)*accum
	return covar

print calc_covar(0,0.5,1.0,qsp,kca)


radii=np.asarray([0.0,1.0,2.0,3.0,4.0])
nradii=len(radii)

nmodes=3

cov_dict=dict()

for l in range(0,nmodes):
	cov_dict[l] = np.zeros((nradii,nradii))
	for ri1 in range(nradii):
		(cov_dict[l])[ri1,ri1] = calc_var(l,radii[ri1],qsp,hca)
		for ri2 in range(ri1+1,nradii):
			(cov_dict[l])[ri1,ri2] = calc_covar(l,radii[ri1],radii[ri2],qsp,kca)
			(cov_dict[l])[ri2,ri1] = (cov_dict[l])[ri1,ri2]

