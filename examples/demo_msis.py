#!/usr/bin/env python3
"""
this program demos the customized MSIS used by GLOW
"""

from dateutil.parser import parse
from numpy import arange,array
#
from glowaurora.runmsisGLOW import rungtdGLOW
#
tselecopts = array([1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],float)

def test_glowmsis(dtime,altkm,glat,glon,f107a,f107,ap,mass):

    glowdens,glowtemp = rungtdGLOW(dtime,altkm,glat,glon,f107a,f107,ap,mass,tselecopts)
#%% now use msise00
    from msise00.runmsis import rungtd1d
    dens,temp = rungtd1d(dtime,altkm,glat,glon,p.f107a,p.f107,p.ap,p.mass,tselecopts)

    #FIXME dim labels the same
    assert (dens.values == glowdens.values).all().all()
    assert (temp.values == glowtemp.values).all().all()
    print('OK')

    return glowdens,glowtemp

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='calls MSISE-00 from Python, a basic demo')
    p.add_argument('simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',type=str,nargs='?',default='2013-04-14T00:00:00Z')
    p.add_argument('-a','--altkm',help='altitude (km) (start,stop,step)',type=float,nargs='+',default=[None])
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(70,0))
    p.add_argument('--f107a',help=' 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)',type=float,default=150)
    p.add_argument('--f107',help='DAILY F10.7 FLUX FOR PREVIOUS DAY',type=float,default=150)
    p.add_argument('--ap',help='daily ap, 0-3hr, 3-6hr, 6-9hr, 9-12hr,12-33hr, 36-57hr',type=float,nargs=7,default=[4,4,4,4,4,4,4])
    p.add_argument('--mass',help=('MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS ' +
                       'CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL. '+
                         'MASS 17 IS Anomalous O ONLY.'),type=float,default=48)
    p = p.parse_args()

    dtime= parse(p.simtime)
#%% altitude 1-D mode
    if p.altkm[0] is None:
        amm = (60,1000,5)
    elif len(p.altkm) == 3:
        amm = p.altkm[0],p.altkm[1],p.altkm[2]
    altkm = arange(amm[0],amm[1],amm[2])
    glat,glon=p.latlon

    print('using altitudes from {:.1f} to {:.1f} km'.format(altkm[0],altkm[-1]))

    dens,temp = test_glowmsis(dtime,altkm,glat,glon,p.f107a,p.f107,p.ap,p.mass)