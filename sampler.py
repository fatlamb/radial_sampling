#! /usr/bin/python

import numpy as np
from numpy.random import multivariate_normal
from numpy import linalg as LA

f = np.load("covariance.npz")

radii = f["radii"]
nradii=len(radii)
#sample = multivariate_normal(f["reduced "

sample_dict=dict()

#Sample the 00 field.
mean = f["reduced_mean"][0]
cov = f["reduced_cov"][0]


w,v = LA.eig(f["cov_dict"][()][4])
print w
print v
#print cov
part1=np.sqrt(4.0*np.pi)*f["nu"]
part2=multivariate_normal(mean,cov)
total = np.zeros(nradii)
total[0] = part1
total[1:] = part2
sample_dict[(0,0)] = total 

#Sample the 1m field.
mean = f["reduced_mean"][1]
cov = f["reduced_cov"][1]
for m in range(-1,2):
	part1=0.0
	part2=multivariate_normal(mean,cov)
	total = np.zeros(nradii)
	total[0] = part1
	total[1:] = part2
	sample_dict[(1,m)] = total 

for l in range(2,f['nmodes']):
	for m in range(-l,l+1):
		sample_dict[(l,m)] = multivariate_normal(np.zeros(nradii),f["cov_dict"][()][l])

