#!/usr/bin/env python3
"""
http://www.netlib.org/toms/493.gz
before running this, type in Terminal (one-time):
f2py -m quartic -h quartic.pyf 493.f
f2py -c quartic.pyf 493.f
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
from numpy import zeros,array
import sys
sys.path.append('fortrancode')
#
from quartic import rpoly

Lin=101

def runquartic(pin):

    polyin = zeros(Lin)
    polyin[:pin.size] = pin
    return rpoly(polyin,pin.size)





if __name__ == '__main__':
    r,i,f = runquartic(array([1,2,3,4]))