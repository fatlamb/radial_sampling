#!/usr/bin/python

from scipy.interpolate import PPoly, splrep
from Crawlers import KCrawlerDict
from Crawlers import HCrawlerDict
from Crawlers import XYCrawlerDict
from scipy.special import spherical_jn
from scipy import integrate


import numpy as np

k,p = np.loadtxt('spectrum.csv',dtype=(float,float), delimiter=",", unpack=True)

#Fit a quadratic spline to the data.
qs = splrep(k,p,k=2,s=0.0)

#Convert this spline to a piecewise quadratic
qsp= PPoly.from_spline(qs,extrapolate=True)

def AboveFirstRoot(x,l):
	#Using a linear interpolation of j_l(x) roots for l between 0 and 100
	#to determine whether the j_l(x) argument lies above the first root.
	if (x > (4.75 + 1.05*l)):
		return True
	else:
		return False

def HZero(x,l,n):
	if l!=0:
		return 0.0
	else:
		return (1.0/(n+1))*x**(n+1)

#Generate a dictionary of XYn. This will fill as Xn/Yn(x) are requested.
#These calls are also idempotent, so no values are calculated twice.	
xyca = XYCrawlerDict()

#Generate a dictionary of Kln. 
hca = HCrawlerDict(xyca)

def HInt(x,l,n):
	return x**n*spherical_jn(l,x)**2

def calc_var(l,r_alpha,local_qsp,local_ca,rzero=False):

	#Extract breakpoints from the piecewise polynomial.
	bkpts = local_qsp.x
	
	#Extract polynomial coefficients.
	coeffs = local_qsp.c
	
	npoints = len(bkpts)
	delta_h = np.zeros((3,npoints-1))
		
#optimization FIXME
	for i in range(npoints-1):
		for n in range(3):
			argpair = [r_alpha*bkpts[i],r_alpha*bkpts[i+1]]
			intvals=[0.0,0.0]
			if (AboveFirstRoot(argpair[0],l)==True and AboveFirstRoot(argpair[1],l)==True):
				for j in range(2):
					if rzero:
						intvals[j] = HZero(argpair[j],l,n)
					else:
						intvals[j]= (1/r_alpha)**(n+1)* local_ca.GetEntry((argpair[j]),(l,n))	
				delta_h[n,i] = intvals[1] - intvals[0]
			else:
				res, err = integrate.quad(HInt, argpair[0], argpair[1], args=(l,n))
				delta_h[n,i]= (1.0/r_alpha)**(n+1.0)*res
	accum = 0
	for i in range(npoints-1):
		accum += (coeffs[0,i] - coeffs[1,i]*bkpts[i] + coeffs[2,i]*bkpts[i]**2)*delta_h[0,i]
		accum += (coeffs[1,i] - 2*coeffs[2,i]*bkpts[i]**2)*delta_h[1,i]
		accum += coeffs[2,i]*delta_h[2,i]
	
	var = 16.0*np.pi**2/(r_alpha**l * r_alpha**l)*accum
	return var



#Generate a dictionary of Kln. 
kca = KCrawlerDict(xyca)

def KInt(x,l,n,alpha,beta):
	return x**n*spherical_jn(l,x*alpha)*spherical_jn(l,x*beta)


def calc_covar(l,r_alpha,r_beta,local_qsp,local_ca):

	#Extract breakpoints from the piecewise polynomial.
	bkpts = local_qsp.x
	
	#Extract polynomial coefficients.
	coeffs = local_qsp.c
	
	npoints = len(bkpts)
	delta_k = np.zeros((3,npoints-1))
	
#optimization FIXME
	for i in range(npoints-1):
		for n in range(3):
			argpair = [bkpts[i],bkpts[i+1]]
			intvals=[0.0,0.0]
			if (AboveFirstRoot(argpair[0],l)==True and AboveFirstRoot(argpair[1],l)==True):
				for j in range(2):
					intvals[j]= local_ca.GetEntry((bkpts[j],r_alpha,r_beta),(l,n))	
				delta_k[n,i] = intvals[1] - intvals[0]
			else:
				delta_k[n,i], err = integrate.quad(KInt, argpair[0], argpair[1], args=(l,n,r_alpha,r_beta))
	
	accum = 0
	for i in range(npoints-1):
		accum += (coeffs[0,i] - coeffs[1,i]*bkpts[i] + coeffs[2,i]*bkpts[i]**2)*delta_k[0,i]
		accum += (coeffs[1,i] - 2*coeffs[2,i]*bkpts[i]**2)*delta_k[1,i]
		accum += coeffs[2,i]*delta_k[2,i]
	
	covar = 16.0*np.pi**2/(r_alpha**l * r_beta**l)*accum
	return covar

print calc_covar(0,0.5,1.0,qsp,kca)


radii=np.asarray([0.0,1.0,2.0,3.0,4.0])
#radii=np.asarray([1.0,2.0,3.0,4.0])
nradii=len(radii)

nmodes=3

cov_dict=dict()

for l in range(0,nmodes):
	cov_dict[l] = np.zeros((nradii,nradii))
	for ri1 in range(nradii):
		if ri1==0:
			(cov_dict[l])[ri1,ri1] = calc_var(l,radii[ri1],qsp,hca,rzero=True)
		else:
			(cov_dict[l])[ri1,ri1] = calc_var(l,radii[ri1],qsp,hca)
		for ri2 in range(ri1+1,nradii):
			(cov_dict[l])[ri1,ri2] = calc_covar(l,radii[ri1],radii[ri2],qsp,kca)
			(cov_dict[l])[ri2,ri1] = (cov_dict[l])[ri1,ri2]

