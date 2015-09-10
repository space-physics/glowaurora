#!/usr/bin/env python3
"""
Registration testing of GLOW
Michael Hirsch

f2py -m glowfort -c egrid.f maxt.f glow.f vquart.f gchem.f ephoto.f solzen.f rcolum.f etrans.f exsect.f ssflux.f snoem.f snoemint.f geomag.f nrlmsise00.f qback.f fieldm.f iri90.f aurora_sub.f --quiet

"""
from datetime import datetime
from numpy import array,zeros,float32,log,arange,append,isclose
from numpy.testing import assert_allclose
from os import chdir,environ
#
from histutils.fortrandates import datetime2yd,datetime2gtd
from msise00.runmsis import rungtd1d
#################################
#TODO hack for module data path issue
chdir(environ['HOME'])
import glowaurora
from glowaurora import glowfort
chdir(glowaurora.__path__[0])
#################################
#%% test inputs
z = arange(80,110+1,1)
z = append(z,array([111.5,113.,114.5,116.,118.,120.,122.,124.,126., 128.,130.,132.,134.,136.,138.,140.,142.,144.,146., 148.,150.,153.,156.,159.,162.,165.,168.,172.,176., 180.,185.,190.,195.,200.,205.,211.,217.,223.,230.,237.,244.,252.,260.,268.,276.,284.,292.,300.,309., 318.,327.,336.,345.,355.,365.,375.,385.,395.,406., 417.,428.,440.,453.,467.,482.,498.,515.,533.,551., 570.,590.,610.,630.,650.,670.,690.,710.,730.,750., 770.,790.,810.,830.,850.,870.,890.,910.,930.,950.]))
nbins = 190; jmax=120 #glow.h
eflux = 1.
e0 = 1e3
maxind = 112
glat = 65; glon=-148
ap=4; f107=100; f107a=100
nmaj=3; nst=6
dtime = datetime(2013,4,14,8,54,0)
#
yd,utsec = datetime2yd(dtime)[:2]

def test_egrid_maxt():
    ener,dE = glowfort.egrid()
    assert_allclose(ener[[maxind,maxind+10,-1]],[1017.7124,1677.9241,47825.418])
#%% test of maxt
    phi = glowfort.maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)
    assert phi.argmax() == maxind
    assert_allclose(phi[[maxind,maxind+10]],[ 114810.6,97814.438])

#%% test vquart (quartic root) KNOWN DEFECTIVE FORTRAN ALGORITHM
#Aquart = tile([-1,0,0,0,1],(jmax,1))
#qroot = glowfort.vquart(Aquart,1)
#assert_allclose(qroot[0],roots(Aquart[0,-1]))
#Aquart = array([[-1,0,0,0,1],
#                [-1,0,0,1,1]])
#nq = Aquart.shape[0]
#Aquart = tile(Aquart,(jmax//nq,1))
#qroot = glowfort.vquartmod.vquart(Aquart, nq)
#try:
#    assert_allclose(qroot[:nq],
#                    [1,0.8191725133961643])
#except AssertionError as e:
#    print('this mismatch is in discussion with S. Solomon.   {}'.format(e))

def test_solzen():
    sza = glowfort.solzen(yd,utsec,glat,glon)
    assert isclose(sza, 104.68412017822266)
    return sza

def test_snoem():
    doy = datetime2gtd(dtime)[0]
    zno,maglat,nozm = glowfort.snoem(doy,1.75*log(0.4*ap),f107)
    assert_allclose((nozm[12,15],nozm[-2,-1]),(33547142.,  44171752.))
    return nozm

def test_snoemint():
    densd,tempd = rungtd1d(dtime,z,glat,glon,f107a,f107,[ap]*7)
# (nighttime background ionization)
    znoint = glowfort.snoemint(dtime.strftime('%Y%j'),glat,glon,f107,ap,z,tempd['heretemp'])
    assert_allclose(znoint[[28,63]], (1.40939280e+08,   4.21317850e+06))
    return znoint

def test_fieldm():
    xdip,ydip,zdip,totfield,dipang,decl,smodip = glowfort.fieldm(glat,glon%360,z[50])
    assert isclose(xdip,0.10698765516281128)
    assert isclose(totfield,0.532055139541626)
    assert isclose(dipang,76.86974334716797)

def test_ssflux():
    iscale=1; hlybr=0.; hlya=0.; fexvir=0.; heiew=0.; xuvfac=3.
    wave1,wave2,sflux = glowfort.ssflux(iscale,f107,f107a,hlybr,fexvir,hlya,heiew,xuvfac)
    assert_allclose(sflux[[11,23]],(4.27225743e+11,   5.54400400e+07))

def test_rcolum_qback():
    densd,tempd = rungtd1d(dtime,z,glat,glon,f107a,f107,[ap]*7)

    """ VCD: Vertical Column Density """
    sza = test_solzen()
    zcol,zvcd = glowfort.rcolum(sza,z*1e5,densd[['O','O2','N2']].values.T,tempd['heretemp'])
# FIXME these tests were numerically unstable (near infinity values)
    assert isclose(zcol[0,0], 1e30) #see rcolum comments for sun below horizon 1e30
    assert isclose(zvcd[2,5],8.04e+25,rtol=1e-2) #TODO changes a bit between python 2 / 3
#%% skipping EPHOTO since we care about night time more for now
    znoint = test_snoemint()
    # zeros because nighttime
    photoi = zeros((nst,nmaj,jmax),dtype=float32,order='F')
    phono = zeros((nst,jmax),dtype=float32,order='F')
    glowfort.qback(zmaj=densd[['O','O2','N2']].values.T,
                                zno=znoint,
                                zvcd=zvcd,
                                photoi=photoi,phono=phono)
    assert isclose(photoi[0,0,117],1.4616582)
    assert isclose(phono[0,73],0.00030911868)

def test_glow():
    # electron precipitation
    """ First enact "glow" subroutine, which calls QBACK, ETRANS and GCHEM among others
    """
    glowfort.glow() #no args

    #%% ver and constituants
    """
    currently using common block CGLOW, in future use module
    """
    zceta = glowfort.cglow.zceta.T
    zeta = glowfort.cglow.zeta.T[:,:11]
    zcsum = zceta.sum(axis=-1)[:,:11]
    assert_allclose(zcsum,zeta,rtol=1e-6)

if __name__ == '__main__':
    test_solzen()
    test_egrid_maxt()
    test_snoem()
    test_snoemint()
    test_fieldm()
    test_ssflux()
    test_rcolum_qback()
    test_glow()
