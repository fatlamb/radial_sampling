#! /usr/bin/python

import numpy as np
from numpy.random import multivariate_normal

f = np.load("covariance.npz")

nradii=len(f["radii"])

