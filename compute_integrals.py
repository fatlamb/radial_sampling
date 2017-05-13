#!/usr/bin/python

import pickle
from scipy.interpolate import CubicSpline
from Crawlers import KCalc
from Crawlers import HCalc
from Crawlers import ICalc
import numpy as np
from tqdm import tqdm
from timeit import default_timer as timer
import os

intpath=os.environ['INTDIR']

spline_degree=3 #Cubic splines

def calc_diag_delta(l,r_alpha,bkpts):
	hcalc = HCalc.HCalc()
	icalc = ICalc.ICalc()
	npoints = len(bkpts)
	delta_h = np.zeros((spline_degree+1,npoints-1))

	for i in range(npoints-1):
		for n in range(spline_degree+1):
			kpair = [bkpts[i],bkpts[i+1]]
			delta_h[n,i] = hcalc.Calculate(kpair,l,n+2,r_alpha)

	return delta_h


def calc_offdiag_delta(l,r_alpha,r_beta,bkpts):
	kcalc = KCalc.KCalc()
	npoints = len(bkpts)
	delta_k = np.zeros((spline_degree+1,npoints-1))

	for i in range(npoints-1):
		for n in range(spline_degree+1):
			kpair = [bkpts[i],bkpts[i+1]]
			delta_k[n,i] = kcalc.Calculate(kpair,l,n+2,r_alpha,r_beta)	
	return delta_k


breakpoints=np.logspace(-4,3,50)
radii=np.linspace(0,5,200)
nradii=len(radii)

modes=range(0,3)

start=timer()
times=[]

for l in tqdm(modes):
	delta_dict=dict()

	for ri1 in range(nradii):
		delta_dict[(ri1,ri1)] = calc_diag_delta(l,radii[ri1],breakpoints)
		for ri2 in range(ri1+1,nradii):
			delta_dict[(ri1,ri2)] = calc_offdiag_delta(l,radii[ri1],radii[ri2],breakpoints)

#	np.savez("intpath+"/integral_parameters"+str(l),bkpts=breakpoints, l=l, radii=radii)
#	pickle.dump(delta_dict, open(intpath+"/integrals"+str(l)+".p","wb"))

	end=timer()
	times.append(end-start)
	start=end

print times
print radii

