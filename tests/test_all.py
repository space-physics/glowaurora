#!/usr/bin/env python
"""
Registration testing of GLOW
Michael Hirsch

f2py -m glowfort -c egrid.f maxt.f glow.f vquart.f gchem.f ephoto.f solzen.f
rcolum.f etrans.f exsect.f ssflux.f snoem.f snoemint.f geomag.f nrlmsise00.f
 qback.f fieldm.f iri90.f aurora_sub.f --quiet

pytest -v
"""
import pytest
from pytest import approx
from datetime import datetime
import numpy as np
import glowfort

# %% test inputs
z = np.arange(30, 110 + 1, 1)
z = np.append(z, [111.5, 113., 114.5, 116.])
z = np.append(z, np.arange(118, 150 + 2, 2))
z = np.append(z, np.arange(153, 168 + 3, 3))
z = np.append(z, np.arange(172, 180 + 4, 4))
z = np.append(z, np.arange(185, 205 + 5, 5))
z = np.append(z, np.arange(211, 223 + 6, 6))
z = np.append(z, np.arange(230, 244 + 7, 7))
z = np.append(z, np.arange(252, 300 + 8, 8))
z = np.append(z, np.arange(309, 345 + 9, 9))
z = np.append(z, np.arange(355, 395 + 10, 10))
z = np.append(z, np.arange(406, 428 + 11, 11))
z = np.append(z, [440, 453, 467, 482, 498, 515, 533, 551])
z = np.append(z, np.arange(570, 950 + 20, 20))

nbins = 190

eflux = 1.
e0 = 1e3
maxind = 112
glat = 70
glon = 0  # like aurora.in
ap = 4
f107 = 100
f107a = 100
dtime = datetime(1999, 12, 21)


def test_egrid_maxt():
    ener, dE = glowfort.egrid()
    assert ener[[maxind, maxind + 10, -1]] == approx([1017.7124, 1677.9241, 47825.418], rel=0.001)
# %% test of maxt
    phi = glowfort.maxt(eflux, e0, ener, dE, itail=0, fmono=np.nan, emono=np.nan)
    assert phi.argmax() == maxind
    assert phi[[maxind, maxind + 10]] == approx([114810.6, 97814.438])

# %% test vquart (quartic root) KNOWN DEFECTIVE FORTRAN ALGORITHM
# Aquart = tile([-1,0,0,0,1],(jmax,1))
# qroot = glowfort.vquart(Aquart,1)
# assert qroot[0] == approx(roots(Aquart[0,-1]))
# Aquart = array([[-1,0,0,0,1],
#                [-1,0,0,1,1]])
# nq = Aquart.shape[0]
# Aquart = tile(Aquart,(jmax//nq,1))
# qroot = glowfort.vquartmod.vquart(Aquart, nq)
# try:
#    assert qroot[:nq] == approx([1,0.8191725133961643])
# except AssertionError as e:
#    print('this mismatch is in discussion with S. Solomon.   {}'.format(e))


def test_fieldm():
    xdip, ydip, zdip, totfield, dipang, decl, smodip = glowfort.fieldm(dlat=glat, dlong=glon % 360, alt=z[50])
    assert xdip == approx(0.1049523800611496)
    assert totfield == approx(0.5043528079986572)
    assert dipang == approx(77.72911071777344)


def test_ssflux():
    iscale = 1
    hlybr = 0.
    hlya = 0.
    fexvir = 0.
    heiew = 0.
    xuvfac = 3.
    wave1, wave2, sflux = glowfort.ssflux(iscale, f107, f107a, hlybr, fexvir, hlya, heiew, xuvfac)
    assert sflux[[11, 23]] == approx((4.27225743e+11, 5.54400400e+07))


# def test_eigen():
#    ener,dE = glowfort.egrid()
#    sim = makeeigen(ener,ones_like(ener),dtime,(glat,glon))


if __name__ == '__main__':
    pytest.main(['-xv', __file__])
