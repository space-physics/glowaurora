#!/usr/bin/env python
import glowaurora as ga
from pytest import approx
import pytest
import sciencedates as sd
from datetime import datetime
import numpy as np

dtime = datetime(1999, 12, 21)
ap = 4
glat = 70
glon = 0  # like aurora.in
f107 = f107a = 100
nmaj = 3
nst = 6
jmax = 170  # glow.h
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

yd, utsec = sd.datetime2yeardoy(dtime)[:2]


def solzen():
    sza = ga.glowfort.solzen(yd, utsec, glat, glon)
    assert sza == approx(133.43113708496094)

    return sza


def snoem():
    doy = sd.datetime2gtd(dtime)[0]

    zno, maglat, nozm = ga.glowfort.snoem(doy, 1.75 * np.log(0.4 * ap), f107)
    assert (nozm[12, 15], nozm[-2, -1]) == approx(35077728.0, 1.118755e+08)

    return nozm


def snoemint():
    msise00 = pytest.importorskip('msise00')

    atmos = msise00.rungtd1d(dtime, z, glat, glon)
# (nighttime background ionization)
    znoint = ga.glowfort.snoemint(dtime.strftime('%Y%j'), glat, glon, f107, ap, z,
                                  atmos['Tn'])
    assert znoint[[28, 143]] == approx((1.262170e+08, 1110.28), rel=1e-5)  # arbitrary
    return znoint


def test_rcolum_qback():
    msise00 = pytest.importorskip('msise00')

    atmos = msise00.rungtd1d(dtime, z, glat, glon)

    """ VCD: Vertical Column Density """
    sza = solzen()
    dens = atmos[['O', 'O2', 'N2']].to_array().values.squeeze()
    zcol, zvcd = ga.glowfort.rcolum(sza, z * 1e5,
                                    dens,
                                    atmos['Tn'])
# FIXME these tests were numerically unstable (near infinity values)
    # see rcolum comments for sun below horizon 1e30
    assert zcol[0, 0] == approx(1e30)
    # TODO changes a bit between python 2 / 3
    assert zvcd[2, 5] == approx(5.97157e+28, rel=1e-2)
# %% skipping EPHOTO since we care about night time more for now
    znoint = snoemint()
    # zeros because nighttime
    photoi = np.zeros((nst, nmaj, jmax), dtype=np.float32, order='F')
    phono = np.zeros((nst, jmax), dtype=np.float32, order='F')

    ga.glowfort.qback(zmaj=dens,
                      zno=znoint,
                      zvcd=zvcd,
                      photoi=photoi, phono=phono)
    # arbitrary point check
    assert photoi[0, 0, 77] == approx(1.38091e-18, rel=1e-5)
    assert phono[0, 73] == approx(0.0, rel=1e-5)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
