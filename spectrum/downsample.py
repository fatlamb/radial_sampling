#!/usr/bin/python

import numpy as np
import scipy.interpolate
from scipy.interpolate import CubicSpline

k,p = np.loadtxt('outspec.csv',dtype=(float,float), delimiter=",", unpack=True)
csp=CubicSpline(k,p)

nsamp=50
kdown = np.logspace(-4,3,nsamp)
pdown = csp(kdown)

deriv= (csp.derivative(nu=1))(k)
#Setting boundary derivatives to match the finer spline.
#For some reason, the default boundary derivatives caused trouble.
cspdown = CubicSpline(kdown,pdown,bc_type=((1,deriv[0]),(1,deriv[-1])))

np.savez("spline_down", bkpts = cspdown.x, coeffs = cspdown.c)
