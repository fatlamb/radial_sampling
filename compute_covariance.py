#!/usr/bin/python

import pickle
from scipy.interpolate import CubicSpline
from numpy import linalg as LA
from math import factorial
from Crawlers import ICalc
import numpy as np
from tqdm import tqdm
from timeit import default_timer as timer
import os

intpath = os.environ['INTDIR']

def calc_var(l,r_alpha,_bkpts,_coeffs,delta_h):

	accum = 0
	for i in range(len(_bkpts)-1):
		accum += (_coeffs[3,i] - _coeffs[2,i]*_bkpts[i] + _coeffs[1,i]*_bkpts[i]**2 - _coeffs[0,i]*_bkpts[i]**3)*delta_h[0,i]
		accum += (_coeffs[2,i] - 2*_coeffs[1,i]*_bkpts[i] + 3*_coeffs[0,i]*_bkpts[i]**2)*delta_h[1,i]
		accum += (_coeffs[1,i] - 3*_coeffs[0,i]*_bkpts[i])*delta_h[2,i]
		accum += (_coeffs[0,i])*delta_h[3,i]

	var = (16.0*np.pi**2)*accum

	if var<0:
		print "NEGAUTOVAR"

	return var

def calc_covar(l,r_alpha,r_beta,_bkpts,_coeffs,delta_k):

	accum = 0

	for i in range(len(_bkpts)-1):
		accum += (_coeffs[3,i] - _coeffs[2,i]*_bkpts[i] + _coeffs[1,i]*_bkpts[i]**2 - _coeffs[0,i]*_bkpts[i]**3)*delta_k[0,i]
		accum += (_coeffs[2,i] - 2*_coeffs[1,i]*_bkpts[i] + 3*_coeffs[0,i]*_bkpts[i]**2)*delta_k[1,i]
		accum += (_coeffs[1,i] - 3*_coeffs[0,i]*_bkpts[i])*delta_k[2,i]
		accum += (_coeffs[0,i])*delta_k[3,i]
	
		covar = (16.0*np.pi**2)*accum

	return covar


spline_data = np.load("spectrum/spline_down.npz")
bkpts = spline_data['bkpts']
coeffs = spline_data['coeffs']

integral_parameters = np.load(intpath+"/integral_parameters0.npz")
#modes = integral_parameters['modes']
radii = integral_parameters['radii']
nradii=len(radii)

modes=range(0,10)
print "MODES: ",modes




cov_dict=dict()

start=timer()
times=[]
#for l in tqdm(range(0,nmodes)):
for l in tqdm(modes):
	delta_dict = pickle.load(open(intpath+"/integrals"+str(l)+".p","rb"))
	cov_dict[l] = np.zeros((nradii,nradii))
	for ri1 in range(nradii):
		(cov_dict[l])[ri1,ri1] = calc_var(l,radii[ri1],bkpts,coeffs,delta_dict[(ri1,ri1)])
		for ri2 in range(ri1+1,nradii):
			(cov_dict[l])[ri1,ri2] = calc_covar(l,radii[ri1],radii[ri2],bkpts,coeffs,delta_dict[(ri1,ri2)])
			(cov_dict[l])[ri2,ri1] = (cov_dict[l])[ri1,ri2]
	end=timer()
	times.append(end-start)
	start=end

print "Completion time[mode]: ",times
#print radii
#print cov_dict

np.savez("covariances/hires_covariance",modes=modes, radii=radii)
pickle.dump(cov_dict, open("covariances/hires_covariance.p","wb"))
