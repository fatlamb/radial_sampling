#! /usr/bin/python
import matplotlib.pyplot as plt
import pickle
import numpy as np
import math
from numpy.random import multivariate_normal
from numpy import linalg as LA
import csv

from pylab import rcParams
rcParams['figure.figsize'] = 16, 8
plt.rc('font',family='Times New Roman')
plt.rc('font', size=16)

def readcsv(fname, curve):
	with open(fname, 'rb') as curvefile:
		curvereader = csv.reader(curvefile, delimiter=',')
		for row in curvereader:
			numarray = np.zeros(len(row))
			for i,el in enumerate(row):
				numarray[i]=el
			curve.append(numarray)

meancurve = list()
readcsv('meangood.csv',meancurve)
locurve = list()
readcsv('loboundgood.csv',locurve)
hicurve = list()
readcsv('hiboundgood.csv',hicurve)




cov_dict = pickle.load(open("hires_covariance.p","rb"))
f = np.load("hires_covariance.npz")

deriv_cov_dict = pickle.load(open("hires_deriv_covariance.p","rb"))
df = np.load("hires_deriv_covariance.npz")

radii = f["radii"]
nradii=len(radii)
#sample = multivariate_normal(f["reduced "
#cov_dict = f["cov_dict"]

#print cov_dict 

sample_dict=dict()
phi_list=[]


#Scrub negative eigenvalues
cov_scrubbed=dict()

def scrubcov(_cov):
	w,v = LA.eig(_cov)
	for n,el in enumerate(w):
		if el<0:
			w[n]=0
	return np.dot(v,np.diag(np.sqrt(w)))

for l in range(0,len(f['modes'])):
	cov_scrubbed[l] = scrubcov(cov_dict[l])


nu=20
bias_val=np.sqrt(4.0*np.pi)*nu
reduced_cov = list()
reduced_cov_scrubbed = list()
reduced_mean = list()
s11 = cov_dict[0][0,0]
s12 = cov_dict[0][0,1:]
s21 = cov_dict[0][1:,0]
s22 = cov_dict[0][1:,1:]

reduced_mean.append(s21*(1.0/s11)*bias_val)
reduced_cov.append(s22 - np.outer(s12,s21)*(1.0/s11))


#Gradiant Biasing
s11 = deriv_cov_dict[-1]
s12 = deriv_cov_dict[1]
s21 = deriv_cov_dict[1]
s22 = cov_dict[1]

reduced_mean.append(np.zeros(nradii))
reduced_cov.append(s22 - np.outer(s12,s21)*(1.0/s11))



for i in range(0,2):
	reduced_cov_scrubbed.append(scrubcov(reduced_cov[i]))


for l in f['modes']:
	negs=[]
	#print "Mode=",l
	w,v = LA.eig(cov_dict[l])
	#print w
	for x in w:
		if x<0:
			negs.append(x)
	#print "ABS: ",np.max(np.abs(negs))
#print v


print "RAD: ", radii

print "DDVAR:",deriv_cov_dict[-1]
print "DCOV1:",deriv_cov_dict[1]

print "COV0", cov_dict[0]
print "COV8", cov_dict[8]
print "COV9", cov_dict[9]


nplots=1
for j in range(0,nplots):
	#Sample the 00 field.
	#print "bleat1"
	mean = reduced_mean[0]
	cov = reduced_cov_scrubbed[0]
	#print mean.shape
	#print cov.shape
		
	z = multivariate_normal(np.zeros(len(mean)),np.diag(np.ones(len(mean))))
	#print "MEAN: ",mean
	#print "Z: ", z
	part2 = mean + np.dot(cov,z)
	part1=bias_val
	#print "P2: ",part2
		#	part2=multivariate_normal(mean,cov)
	#print "bleat1"
	total = np.zeros(nradii)
	total[0] = part1
	total[1:] = part2
	sample_dict[(0,0)] = total 
	#print "Part1: ",part1
	#print "TOTAL: ", total
	
	#Sample the 1m field.
	#print "bleat2"
	mean = reduced_mean[1]
	cov = reduced_cov_scrubbed[1]
	for m in range(-1,2):
		z = multivariate_normal(np.zeros(len(mean)),np.diag(np.ones(len(mean))))
		sample_dict[(1,m)] = mean + np.dot(cov,z)
		#part2 = mean + np.dot(cov,z)
		#part1=0.0
		##part2=multivariate_normal(mean,cov)
		#total = np.zeros(nradii)
		#total[0] = part1
		#total[1:] = part2
		#sample_dict[(1,m)] = total 
	#print "bleat3"
	#for l in range(0,len(f['modes'])):
	for l in range(2,len(f['modes'])):
		cov=cov_scrubbed[l]
		for m in range(-l,l+1):
			z = multivariate_normal(np.zeros(nradii),np.diag(np.ones(nradii)))
			sample_dict[(l,m)] = np.dot(cov,z)
		
	Phi = np.zeros(nradii)
	l_stack=list()
	for l in range(0,len(f['modes'])):
		l_stack.append(np.zeros(nradii))
	for l in range(0,len(f['modes'])):
	#	for l in range(0,2):
		for m in range(-l,l+1):
			for r in range(nradii):
				Phi[r] += sample_dict[(l,m)][r]**2
				l_stack[l][r] = Phi[r]
	Phi/=(4*np.pi)
	l_stack[l]/=(4*np.pi)
	phi_list.append(Phi)
"""
colors=['red','orange','yellow','blue','purple']

for l in range(0,4):
	fig,ax = plt.subplots()
	ax.set_ylabel(r'$\phi_{lm}$')
	ax.set_xlabel(r'$r$')
	#ax.plot(radii,sample_dict[(1,0)])
	label = r'$\ell = $'+str(l+5)
	ax.plot(radii,sample_dict[(l+5,0)],color=colors[l],label=label)
	ax.legend(loc='upper right',prop={'size':10}, shadow=False,fancybox=True)

	plt.savefig("hires_modeplots/modeplot_l="+str(l+5)+".png")
"""

for i in range(0,nplots):
	fig,ax = plt.subplots()
	ax.set_ylabel(r'$\Phi$')
	plt.title(r'$N=1$'+ " "+r'$m_0=5$'+" "+r'$m_\psi$=0.1'+" "+r'$\nu=20$')
	
	ax.set_xlabel(r'$r$')
	ax.plot(locurve[0],locurve[1],lw=2,label=r'$+1\sigma$')
	ax.plot(meancurve[0],meancurve[1],lw=2,label=r'$\langle\Phi\rangle$')
	ax.plot(hicurve[0],hicurve[1],lw=2,label=r'$+1\sigma$')
#	for k in range(0,i):
#		ax.plot(radii,phi_list[k],'k')
	ax.plot(radii,Phi,'k')
	ax.plot(radii[0],Phi[0],'ko',ms=6)
	ax.legend(loc='upper right',prop={'size':16}, shadow=False,fancybox=True)
#	plt.savefig("new_env_stack/env"+str(i)+".png")
	plt.show()


#plt.savefig("hires_modeplots/modestack.png")
#ax.set_yscale('log')
#plt.show()

#np.savez("radial_sample",radii=radii, sample_dict=sample_dict, nu=nu,phi_list=phi_list)
