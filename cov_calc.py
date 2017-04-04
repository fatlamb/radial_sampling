#!/usr/bin/python

from numpy import linalg as LA
from scipy.interpolate import PPoly, splrep
from math import factorial
from Crawlers import KCalc
from Crawlers import HCalc
from Crawlers import ICalc
import numpy as np

k,p = np.loadtxt('spectrum.csv',dtype=(float,float), delimiter=",", unpack=True)

#Fit a quadratic spline to the data.
qs = splrep(k,p,k=2,s=0.0)

#Convert this spline to a piecewise quadratic
qsp= PPoly.from_spline(qs,extrapolate=True)


def calc_var(l,r_alpha,local_qsp,hcalc,icalc):
	#Extract breakpoints from the piecewise polynomial.
	bkpts = local_qsp.x
	#Extract polynomial coefficients.
	coeffs = local_qsp.c
	npoints = len(bkpts)
	delta_h = np.zeros((3,npoints-1))

	zflag = icalc.DetZero(r_alpha)

	#optimization FIXME
	for i in range(npoints-1):
		for n in range(3):
			kpair = [bkpts[i],bkpts[i+1]]

			if zflag:
				delta_h[n,i] = (1.0/(3.0+2*l+n))*(kpair[1]**(3+2*l+n) - kpair[0]**(3+2*l+n))
			if (not zflag):
				delta_h[n,i] = hcalc.Calculate(kpair,l,n+2,r_alpha)


#	for i in range(npoints-1):
#		accum += (coeffs[0,i] - coeffs[1,i]*bkpts[i] + coeffs[2,i]*bkpts[i]**2)*delta_h[0,i]
#		accum += (coeffs[1,i] - 2*coeffs[2,i]*bkpts[i]**2)*delta_h[1,i]
#		accum += coeffs[2,i]*delta_h[2,i]

#Summation is separated for future precalculation of integrals, so we can rapidly compare power spectra.
	accum = 0
	for i in range(npoints-1):
		for n in range(3):
			accum += (coeffs[n,i])*delta_h[n,i]


	if zflag:
		var = (16.0*np.pi**2)*2**(2*l+2)*factorial(l+1)**2/(factorial(2*(l+1)))**2*accum
	if (not zflag):
		var = (16.0*np.pi**2/(r_alpha**l*r_alpha**l))*accum

	if var<0:
		print "NEGAUTOVAR"

	return var


def calc_covar(l,r_alpha,r_beta,local_qsp,kcalc,icalc):
	#Extract breakpoints from the piecewise polynomial.
	bkpts = local_qsp.x
	#Extract polynomial coefficients.
	coeffs = local_qsp.c
	npoints = len(bkpts)
	delta_k = np.zeros((3,npoints-1))

	azflag = icalc.DetZero(r_alpha)
	bzflag = icalc.DetZero(r_beta)
	
#optimization FIXME
	for i in range(npoints-1):
		for n in range(3):
			kpair = [bkpts[i],bkpts[i+1]]
			if azflag and bzflag:
				delta_k[n,i] = (1.0/(3.0+2*l+n))*(kpair[1]**(3+2*l+n) - kpair[0]**(3+2*l+n))
			if azflag and (not bzflag):
				delta_k[n,i] = icalc.Calculate(kpair,l,n+l+2,r_beta)
			if (not azflag) and bzflag:
				delta_k[n,i] = icalc.Calculate(kpair,l,n+l+2,r_alpha)
			if (not azflag) and (not bzflag):
				delta_k[n,i] = kcalc.Calculate(kpair,l,n+2,r_alpha,r_beta)	

#			if(delta_k[n,i]<0):
#				print "KBAD!!!"

	#Summation is separated for future precalculation of integrals, so we can rapidly compare power spectra.
	accum = 0
	for i in range(npoints-1):
		for n in range(3):
			accum += (coeffs[n,i])*delta_k[n,i]

#	for i in range(npoints-1):
#		accum += (coeffs[0,i] - coeffs[1,i]*bkpts[i] + coeffs[2,i]*bkpts[i]**2)*delta_k[0,i]
#		accum += (coeffs[1,i] - 2*coeffs[2,i]*bkpts[i]**2)*delta_k[1,i]
#		accum += coeffs[2,i]*delta_k[2,i]
	
	if azflag and bzflag:
		covar = (16.0*np.pi**2)*2**(2*l+2)*factorial(l+1)**2/(factorial(2*(l+1)))**2*accum
	if azflag and (not bzflag):
		covar = (16.0*np.pi**2/r_beta**l)*2**(l+1)*factorial(l+1)/(factorial(2*(l+1)))*accum
	if (not azflag) and bzflag:
		covar = (16.0*np.pi**2/r_alpha**l)*2**(l+1)*factorial(l+1)/(factorial(2*(l+1)))*accum
	if (not azflag) and (not bzflag):
		covar = (16.0*np.pi**2/(r_alpha**l*r_beta**l))*accum

	return covar

#print calc_covar(0,0.5,1.0,qsp)
#print calc_var(0,0.5,qsp)

radii=np.asarray([0.0,1.0,2.0])
#radii=np.asarray([1.0,2.0,3.0,4.0])
nradii=len(radii)

nmodes=5

cov_dict=dict()
hcalc = HCalc.HCalc()
kcalc = KCalc.KCalc()
icalc = ICalc.ICalc()

for l in range(0,nmodes):
	cov_dict[l] = np.zeros((nradii,nradii))
	for ri1 in range(nradii):
		(cov_dict[l])[ri1,ri1] = calc_var(l,radii[ri1],qsp,hcalc,icalc)
		for ri2 in range(ri1+1,nradii):
			(cov_dict[l])[ri1,ri2] = calc_covar(l,radii[ri1],radii[ri2],qsp,kcalc,icalc)
			(cov_dict[l])[ri2,ri1] = (cov_dict[l])[ri1,ri2]
print cov_dict


w,v = LA.eig(cov_dict[0])
print w
print v



#Biasing at origin to specific value and zero gradient!
nu=3
bias_vals=[np.sqrt(4.0*np.pi)*nu,0.0]
reduced_cov = list()
reduced_mean = list()

for l in range(0,2):
	s11 = cov_dict[l][0,0]
	s12 = cov_dict[l][0,1:]
	s21 = cov_dict[l][1:,0]
	s22 = cov_dict[l][1:,1:]

	reduced_mean.append(s21*(1.0/s11)*bias_vals[l])
	reduced_cov.append(s22 - np.outer(s12,s21)*(1.0/s11))

#	s22 = cov_dict[l][-1,-1]
#	s21 = cov_dict[l][0,0:-1]
#	s12 = cov_dict[l][0:-1,0]
#	s11 = cov_dict[l][0:-1,0:-1]

#	reduced_mean.append(s12*(1.0/s22)*bias_vals[l])
#	reduced_cov.append(s11 - np.outer(s12,s21)*(1.0/s22))


np.savez("covariance",nmodes=nmodes, radii=radii, cov_dict=cov_dict,nu=nu,reduced_cov=reduced_cov,reduced_mean=reduced_mean)
