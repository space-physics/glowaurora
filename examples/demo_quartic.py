#!/usr/bin/env python3
"""
http://www.netlib.org/toms/493.gz
before running this, type in Terminal (one-time):

cd ../fortran
f2py -m quartic -c quartic.pyf 493.f

-----------------------
>>> print quartic.rpoly.__doc__
zeror,zeroi,fail = rpoly(op,degree)

Wrapper for ``rpoly``.

Parameters
----------
op : input rank-1 array('d') with bounds (101)
degree : in/output rank-0 array(int,'i')

Returns
-------
zeror : rank-1 array('d') with bounds (100)
zeroi : rank-1 array('d') with bounds (100)
fail : int


"""
from __future__ import division,print_function
from numpy import zeros,repeat,array
import sys
sys.path.append('../fortran')
#
from glowaurora import glowfort
try:
    from quartic import rpoly
except Exception as e:
    print(e)
    print('in Terminal, type: ')
    print('f2py -m quartic -c 493.f')

Lrpoly = 101 #493.f
Lin=170 #must match compiled cglow.h
nj=20 #arbitrary
nc=5 #per vquart.f

def runvquart(pin):
    polyin = repeat(pin[:,None],nc,1) #this is what vquart.f expects
    return glowfort.vquart(polyin,pin.size)

def runrpoly(pin):
    return rpoly(pin,pin.size)

if __name__ == '__main__':
    pin = array([1,0,0,1])
    polyin = zeros(Lin)
    polyin[:pin.size] = pin
#%% defective glow routine
    root = runvquart(polyin)[:pin.size]
#%% correct 493.f routine
    r,i,f = runrpoly(polyin[:Lrpoly])

    print('real {}  \nimag {}'.format(r[:pin.size],i[:pin.size]))
