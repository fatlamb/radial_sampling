#!/usr/bin/python

import pickle
from scipy.interpolate import CubicSpline
from Crawlers import KCalc
from Crawlers import HCalc
from Crawlers import ICalc
import numpy as np
from tqdm import tqdm
from timeit import default_timer as timer

spline_degree=3 #Cubic splines

def calc_diag_delta(l,r_alpha,bkpts):
	hcalc = HCalc.HCalc()
	icalc = ICalc.ICalc()
	npoints = len(bkpts)
	delta_h = np.zeros((spline_degree+1,npoints-1))

	zflag = icalc.DetZero(r_alpha)

	for i in range(npoints-1):
		for n in range(spline_degree+1):
			kpair = [bkpts[i],bkpts[i+1]]

			if zflag:
				delta_h[n,i] = (1.0/(3.0+2*l+n))*(kpair[1]**(3+2*l+n) - kpair[0]**(3+2*l+n))
			if (not zflag):
				delta_h[n,i] = hcalc.Calculate(kpair,l,n+2,r_alpha)

	return delta_h


def calc_offdiag_delta(l,r_alpha,r_beta,bkpts):
	kcalc = KCalc.KCalc()
	icalc = ICalc.ICalc()
	npoints = len(bkpts)
	delta_k = np.zeros((spline_degree+1,npoints-1))

	azflag = icalc.DetZero(r_alpha)
	bzflag = icalc.DetZero(r_beta)

	for i in range(npoints-1):
		for n in range(spline_degree+1):
			kpair = [bkpts[i],bkpts[i+1]]
			if azflag and bzflag:
				delta_k[n,i] = (1.0/(3.0+2*l+n))*(kpair[1]**(3+2*l+n) - kpair[0]**(3+2*l+n))
			if azflag and (not bzflag):
				delta_k[n,i] = icalc.Calculate(kpair,l,n+l+2,r_beta)
			if (not azflag) and bzflag:
				delta_k[n,i] = icalc.Calculate(kpair,l,n+l+2,r_alpha)
			if (not azflag) and (not bzflag):
				delta_k[n,i] = kcalc.Calculate(kpair,l,n+2,r_alpha,r_beta)	
	return delta_k


breakpoints=np.logspace(-4,3,200)
radii=np.linspace(0,1000,200)
nradii=len(radii)

modes=range(0,11)

delta_dict=dict()

start=timer()
times=[]

for l in tqdm(modes):
	for ri1 in range(nradii):
		delta_dict[(l,ri1,ri1)] = calc_diag_delta(l,radii[ri1],breakpoints)
		for ri2 in range(ri1+1,nradii):
			delta_dict[(l,ri1,ri2)] = calc_offdiag_delta(l,radii[ri1],radii[ri2],breakpoints)
	end=timer()
	times.append(end-start)
	start=end

print times
print radii

np.savez("integral_parameters_ft",bkpts=breakpoints, modes=modes, radii=radii)
pickle.dump(delta_dict, open("integrals_ft.p","wb"))
