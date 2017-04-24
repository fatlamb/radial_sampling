#!/usr/bin/python

import pickle
from scipy.interpolate import CubicSpline
from numpy import linalg as LA
from math import factorial
from Crawlers import ICalc
import numpy as np
from tqdm import tqdm
from timeit import default_timer as timer

def calc_var(l,r_alpha,spline,delta_h):
	#Extract breakpoints from the piecewise polynomial.
	bkpts = spline.x
	#Extract polynomial coefficients.
	coeffs = spline.c

	icalc = ICalc.ICalc()
	zflag = icalc.DetZero(r_alpha)

	accum = 0
	for i in range(len(bkpts)-1):
		accum += (coeffs[3,i] - coeffs[2,i]*bkpts[i] + coeffs[1,i]*bkpts[i]**2 - coeffs[0,i]*bkpts[i]**3)*delta_h[0,i]
		accum += (coeffs[2,i] - 2*coeffs[1,i]*bkpts[i] + 3*coeffs[0,i]*bkpts[i]**2)*delta_h[1,i]
		accum += (coeffs[1,i] - 3*coeffs[0,i]*bkpts[i])*delta_h[2,i]
		accum += (coeffs[0,i])*delta_h[3,i]

	if zflag:
		var = (16.0*np.pi**2)*2**(2*l+2)*factorial(l+1)**2/(factorial(2*(l+1)))**2*accum
	if (not zflag):
		var = (16.0*np.pi**2/(r_alpha**l*r_alpha**l))*accum

	if var<0:
		print "NEGAUTOVAR"

	return var

def calc_covar(l,r_alpha,r_beta,spline,delta_k):
	#Extract breakpoints from the piecewise polynomial.
	bkpts = spline.x
	#Extract polynomial coefficients.
	coeffs = spline.c

	icalc = ICalc.ICalc()
	azflag = icalc.DetZero(r_alpha)
	bzflag = icalc.DetZero(r_beta)

	accum = 0

	for i in range(len(bkpts)-1):
		accum += (coeffs[3,i] - coeffs[2,i]*bkpts[i] + coeffs[1,i]*bkpts[i]**2 - coeffs[0,i]*bkpts[i]**3)*delta_k[0,i]
		accum += (coeffs[2,i] - 2*coeffs[1,i]*bkpts[i] + 3*coeffs[0,i]*bkpts[i]**2)*delta_k[1,i]
		accum += (coeffs[1,i] - 3*coeffs[0,i]*bkpts[i])*delta_k[2,i]
		accum += (coeffs[0,i])*delta_k[3,i]
	
	if azflag and bzflag:
		covar = (16.0*np.pi**2)*2**(2*l+2)*factorial(l+1)**2/(factorial(2*(l+1)))**2*accum
	if azflag and (not bzflag):
		covar = (16.0*np.pi**2/r_beta**l)*2**(l+1)*factorial(l+1)/(factorial(2*(l+1)))*accum
	if (not azflag) and bzflag:
		covar = (16.0*np.pi**2/r_alpha**l)*2**(l+1)*factorial(l+1)/(factorial(2*(l+1)))*accum
	if (not azflag) and (not bzflag):
		covar = (16.0*np.pi**2/(r_alpha**l*r_beta**l))*accum

	return covar


#k,p = np.loadtxt('spectrum.csv',dtype=(float,float), delimiter=",", unpack=True)

k=np.logspace(-4,3,200)
#p=0.9*(np.exp(-20*k)+ 0.1*np.exp(-5*k))

p=np.exp(-k)

integral_parameters = np.load("integral_parameters.npz")
modes = integral_parameters['modes']
radii = integral_parameters['radii']
nradii=len(radii)

delta_dict = pickle.load(open("integrals.p","rb"))

#We may have to resample in the future to make the k and bkpts match.
#This shouldn't cause trouble because we can compute the spectra
#at very high resolution without too much trouble.

#Fit a cubic spline to the data.
csp = CubicSpline(k,p)

cov_dict=dict()
print "MODES: ",modes

start=timer()
times=[]
#for l in tqdm(range(0,nmodes)):
for l in tqdm(modes):
	print "l", l
	cov_dict[l] = np.zeros((nradii,nradii))
	for ri1 in range(nradii):
		(cov_dict[l])[ri1,ri1] = calc_var(l,radii[ri1],csp,delta_dict[(l,ri1,ri1)])
		for ri2 in range(ri1+1,nradii):
			(cov_dict[l])[ri1,ri2] = calc_covar(l,radii[ri1],radii[ri2],csp,delta_dict[(l,ri1,ri2)])
			(cov_dict[l])[ri2,ri1] = (cov_dict[l])[ri1,ri2]
	end=timer()
	times.append(end-start)
	start=end

print times
print radii
print cov_dict

np.savez("covariance",modes=modes, radii=radii, cov_dict=cov_dict)
