#!/usr/bin/python

import pickle
from scipy.interpolate import CubicSpline
from scipy import integrate as sciint
import scipy.special as special
import numpy as np
from tqdm import tqdm
from timeit import default_timer as timer
import os
import sys

#Load target path for integral dictionaries. 
intpath=os.environ['INTDIR']

#Load working path from environment.
pbhdir = os.environ['PBH']

#Import the levin module
sys.path.append(pbhdir+'/levin_integrals')
import levin

spline_degree=3 #Cubic splines

def k2_analytic(a,b,alpha,beta,l):
	"""
	Substitutes parameters into the analytic result
	for Integrate(k^2j_l(alpha*k)j_l(beta*k),{k,a,b})	
	that was derived in arXiv:1703.06428
	"""
	res_a = (a**2/(alpha**2-beta**2))*(beta*special.spherical_jn(l,alpha*a)
		*special.spherical_jn(l-1,beta*a) - alpha*special.spherical_jn(l-1,alpha*a)
		*special.spherical_jn(l,beta*a))

	T1 = beta*special.spherical_jn(l,alpha*a)*special.spherical_jn(l-1,beta*a)
	T2 = -1.0*alpha*special.spherical_jn(l-1,alpha*a)*special.spherical_jn(l,beta*a)

	sumsize = abs(T1)+abs(T2)
	diffsize = abs(T1+T2)
	if sumsize!=0:
		print "PrecLoss: ",diffsize/sumsize

	res_b = (b**2/(alpha**2-beta**2))*(beta*special.spherical_jn(l,alpha*b)
		*special.spherical_jn(l-1,beta*b) - alpha*special.spherical_jn(l-1,alpha*b)
		*special.spherical_jn(l,beta*b))
	return res_b - res_a


def k0(k):
	return 1.0
def k1(k):
	return k
def k2(k):
	return k**2
def k3(k):
	return k**3

poly_funcs = [k0,k1,k2,k3]

integrator = levin.SphericalBesselIntegrator()

def calc_relerr(err,result):
	if err!=0.0:
		rel_err = abs(err/result)
	else:
		rel_err = 0
	return rel_err

def calc_diag_delta(l,alpha,bkpts):
	npoints = len(bkpts)
	delta = np.zeros((spline_degree+1,npoints-1))
	for i in range(npoints-1):
		for n in range(spline_degree+1):
			delta[n,i] = integrator.HCalc(bkpts[i],bkpts[i+1],alpha,l,poly_funcs[n])
	return delta

def calc_offdiag_delta(l,alpha,beta,bkpts):
	npoints = len(bkpts)
	delta = np.zeros((spline_degree+1,npoints-1))
#	print "                                alpha,beta: ", alpha[0]," , ",beta[0]
	for i in range(npoints-1):
		for n in range(spline_degree+1):
			delta[n,i] = integrator.KCalc(bkpts[i],bkpts[i+1],alpha,beta,l,poly_funcs[n])
			#print("Numerical: ",'{:3E}'.format(delta[n,i]))
#			if l>0 and n==2:
#				analytic = k2_analytic(bkpts[i],bkpts[i+1],alpha[0],beta[0],l)
#				abserr = abs(delta[n,i]-analytic)				
#				relerr = calc_relerr(abserr,analytic)				
#				if relerr>1e-6:
#					print("Deviation from analytic result greater than 1e-6")
#					print("[",bkpts[i],",",bkpts[i+1],"]")
#					print("Analytic: ",analytic)
#					print("AbsDev: ",abserr)
#					print("RelDev: ",relerr)
	return delta


breakpoints=np.logspace(-4,3,50)
radii=np.linspace(0,5,200)
nradii=len(radii)

modes=range(0,3)

start=timer()
times=[]

for l in tqdm(modes):
	delta_dict=dict()

	for ri1 in range(nradii):
		delta_dict[(ri1,ri1)] = calc_diag_delta(l,[radii[ri1],ri1],breakpoints)
		for ri2 in range(ri1+1,nradii):
			delta_dict[(ri1,ri2)] = calc_offdiag_delta(l,[radii[ri1],ri1],[radii[ri2],ri2],breakpoints)

#	np.savez(intpath+"/integral_parameters"+str(l),bkpts=breakpoints, l=l, radii=radii)
#	pickle.dump(delta_dict, open(intpath+"/integrals"+str(l)+".p","wb"))

	end=timer()
	times.append(end-start)
	start=end

print times
