#!/usr/bin/env python3
"""
Registration testing of GLOW
Michael Hirsch
"""
from datetime import datetime
from fortrandates import datetime2gtd
from numpy import array,tile,roots,log,arange,append
from numpy.testing import assert_allclose
import sys
sys.path.append('../msise-00')
from demo_msis import rungtd1d
#%% test inputs
z = arange(80,110+1,1)
z = append(z,array([111.5,113.,114.5,116.,118.,120.,122.,124.,126., 128.,130.,132.,134.,136.,138.,140.,142.,144.,146., 148.,150.,153.,156.,159.,162.,165.,168.,172.,176., 180.,185.,190.,195.,200.,205.,211.,217.,223.,230.,237.,244.,252.,260.,268.,276.,284.,292.,300.,309., 318.,327.,336.,345.,355.,365.,375.,385.,395.,406., 417.,428.,440.,453.,467.,482.,498.,515.,533.,551., 570.,590.,610.,630.,650.,670.,690.,710.,730.,750., 770.,790.,810.,830.,850.,870.,890.,910.,930.,950.]))
nbins = 190; jmax=120 #glow.h
eflux = 1.
e0 = 1e3
maxind = 112
glat = 65; glon=-148
ap=4; f107=100; f107a=100
dtime = datetime(2013,4,14,8,54,0)
#
iyd,utsec = datetime2gtd(dtime)[:2]
#%% test of egrid
from glowgrid import energygrid,maxt
ener,dE = energygrid(nbins)
assert_allclose(ener[[maxind,maxind+10,-1]],[1017.7124,1677.9241,47825.418])
#%% test of maxt
phi = maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)
assert phi.argmax() == maxind
assert_allclose(phi[[maxind,maxind+10]],[ 114810.6,97814.438])
#%% test vquart (quartic root)
from aurora import vquartmod
Aquart = tile([-1,0,0,0,1],(jmax,1))
qroot = vquartmod(Aquart,1)
assert_allclose(qroot[0],roots(Aquart[0,-1]))
#%% test snoem
from aurora import snoemmod
zno,maglat,nozm = snoemmod(iyd,1.75*log(0.4*ap),f107)
assert_allclose((nozm[12,15],nozm[-2,-1]),(33547142.,  44171752.))
#%% test snoemint
from aurora import snoemint
densd,tempd = rungtd1d(dtime,z,glat,glon,f107a,f107,[ap]*7)
znoint = snoemint(dtime.strftime('%Y%j'),glat,glon,f107,ap,z,tempd['heretemp'])
assert_allclose(znoint[[28,63]], (1.40939280e+08,   4.11383025e+06))