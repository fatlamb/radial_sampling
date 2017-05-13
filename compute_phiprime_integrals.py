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

def calc_deriv_deriv(bkpts):
	npoints = len(bkpts)
	delta_self = np.zeros((spline_degree+1,npoints-1))

	for i in range(npoints-1):
		for n in range(spline_degree+1):
			kpair = [bkpts[i],bkpts[i+1]]
			delta_self[n,i] =  (1.0/(5.0+n))*(kpair[1]**(5+n)-kpair[0]**(5+n))

	return delta_self


def calc_deriv_phi(r_alpha,bkpts):
	icalc = ICalc.ICalc()
	npoints = len(bkpts)
	delta_i = np.zeros((spline_degree+1,npoints-1))

	for i in range(npoints-1):
		for n in range(spline_degree+1):
			kpair = [bkpts[i],bkpts[i+1]]
			delta_i[n,i] = icalc.Calculate(kpair,1,n+3,r_alpha)

	return delta_i


breakpoints=np.logspace(-4,3,50)
radii=np.linspace(0,5,200)
nradii=len(radii)

modes=range(1,2)

start=timer()
times=[]

for l in tqdm(modes):
	delta_dict=dict()

	delta_dict[(0,0)] = calc_deriv_deriv(breakpoints)
	for ri in range(0,nradii):
		delta_dict[(0,ri+1)] = calc_deriv_phi(radii[ri],breakpoints)

	np.savez("scaled_output_05_200/deriv_integral_parameters"+str(l),bkpts=breakpoints, l=l, radii=radii)
	pickle.dump(delta_dict, open("scaled_output_05_200/deriv_integrals"+str(l)+".p","wb"))

	end=timer()
	times.append(end-start)
	start=end

print times
print radii

