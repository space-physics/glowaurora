#!/usr/bin/env python3
"""
http://www.netlib.org/toms/493.gz
before running this, type in Terminal (one-time):

f2py3 -m quartic -c fortran/493.f

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
from numpy import zeros, repeat, array, roots
#
import glowfort
import sys
sys.path.append('../fortran')
try:
    from quartic import rpoly
except Exception as e:
    print(e)
    print('in Terminal, type: ')
    print('f2py3 -m quartic -c fortran/493.f')

Lrpoly = 101  # 493.f
Lin = 170  # must match compiled cglow.h
nj = 20  # arbitrary
nc = 5  # per vquart.f


def rootnumpy(pin):
    R = roots(pin)
    return R[R.imag == 0].real


def runvquart(pin):
    polyin = repeat(pin[:, None], nc, 1)  # this is what vquart.f expects
    return glowfort.vquart(polyin, pin.size)


def runrpoly(pin):
    rootreal, rootimag, fail = rpoly(pin, pin.size)
    if not fail:
        return rootreal


if __name__ == '__main__':
    pin = array([1, 0, 0, 1])
    polyin = zeros(Lin)
    polyin[:pin.size] = pin
# %% defective glow routine
    root = runvquart(polyin)[:pin.size]
# %% correct 493.f routine
    r = runrpoly(polyin[:Lrpoly])
# %% numpy
    rootnumpy = rootnumpy(pin)

    print('vquart.f: {}'.format(root))
    print('493.f: {}'.format(r[:pin.size]))
    print('numpy {}'.format(rootnumpy))
