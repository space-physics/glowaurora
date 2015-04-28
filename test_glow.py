#!/usr/bin/env python3
"""
Registration testing of GLOW
Michael Hirsch
"""
from numpy import array,gradient
from numpy.testing import assert_allclose
#%% test inputs
nbins = 190
eflux = 1.
e0 = 1e3
maxind = 112
#%% test of egrid
from glowgrid import energygrid,maxt
ener,dE = energygrid(nbins)
assert_allclose(ener[[maxind,maxind+10,-1]],[1017.7124,1677.9241,47825.418])
#%% test of maxt
phi = maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)
assert phi.argmax() == maxind
assert_allclose(phi[[maxind,maxind+10]],[ 114810.6,97814.438])
