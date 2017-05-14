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

def calc_deriv_deriv(_bkpts,_coeffs,delta_self):
	accum = 0
	for i in range(len(_bkpts)-1):
		accum += (_coeffs[3,i] - _coeffs[2,i]*_bkpts[i] + _coeffs[1,i]*_bkpts[i]**2 - _coeffs[0,i]*_bkpts[i]**3)*delta_self[0,i]
		accum += (_coeffs[2,i] - 2*_coeffs[1,i]*_bkpts[i] + 3*_coeffs[0,i]*_bkpts[i]**2)*delta_self[1,i]
		accum += (_coeffs[1,i] - 3*_coeffs[0,i]*_bkpts[i])*delta_self[2,i]
		accum += (_coeffs[0,i])*delta_self[3,i]

	var = (((4.0*np.pi)**2)/9.0)*accum

	if var<0:
		print "NEGAUTOVAR"

	return var

def calc_deriv_phi(_bkpts,_coeffs,delta_i):

	accum = 0
	for i in range(len(_bkpts)-1):
		accum += (_coeffs[3,i] - _coeffs[2,i]*_bkpts[i] + _coeffs[1,i]*_bkpts[i]**2 - _coeffs[0,i]*_bkpts[i]**3)*delta_i[0,i]
		accum += (_coeffs[2,i] - 2*_coeffs[1,i]*_bkpts[i] + 3*_coeffs[0,i]*_bkpts[i]**2)*delta_i[1,i]
		accum += (_coeffs[1,i] - 3*_coeffs[0,i]*_bkpts[i])*delta_i[2,i]
		accum += (_coeffs[0,i])*delta_i[3,i]
	
	covar = (((4.0*np.pi)**2)/3.0)*accum

	return covar


spline_data = np.load("spectrum/spline_down.npz")
bkpts = spline_data['bkpts']
coeffs = spline_data['coeffs']

integral_parameters = np.load(intpath+"/integral_parameters0.npz")
#modes = integral_parameters['modes']
radii = integral_parameters['radii']
nradii=len(radii)

modes=range(1,2)
#print "MODES: ",modes

cov_dict=dict()

start=timer()
times=[]
#for l in tqdm(range(0,nmodes)):
for l in tqdm(modes):
	delta_dict = pickle.load(open(intpath+"/deriv_integrals"+str(l)+".p","rb"))
	cov_dict[l] = np.zeros(nradii)
	cov_dict[-1] = calc_deriv_deriv(bkpts,coeffs,delta_dict[(0,0)])
	#print "SHEP: ",delta_dict[(0,1)]
	for ri in range(0,nradii):
		(cov_dict[l])[ri] = calc_deriv_phi(bkpts,coeffs,delta_dict[(0,ri+1)])
	end=timer()
	times.append(end-start)
	start=end

print "Completion time[mode]: ",times
#print radii
#print cov_dict

np.savez("covariances/hires_deriv_covariance",modes=modes, radii=radii)
pickle.dump(cov_dict, open("covariances/hires_deriv_covariance.p","wb"))
