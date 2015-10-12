#!/usr/bin/env python3
"""
this program demos the customized MSIS-E00 used by GLOW
"""
from __future__ import division,print_function
import logging
from numpy import repeat,empty,array
from pandas import DataFrame
from dateutil.parser import parse
from numpy import arange, atleast_1d
from os import chdir
#
import glowaurora
from glowaurora.glowfort import gtd7,tselec,csw,meters
glowpath=glowaurora.__path__[0]
#
from histutils.fortrandates import datetime2gtd

tselecopts = array([1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],float)

def rungtd1d(dtime,altkm,glat,glon,f107a,f107,ap,mass):
    chdir(glowpath)
    ap = atleast_1d(ap)
    if ap.size==1: ap = repeat(ap,7)
    species = ['He','O','N2','O2','Ar','Total','H','N','AnomalousO']
    ttypes = ['exotemp','heretemp']

    tselec(tselecopts) #like the msis_driver example
    logging.debug('tselec options used:   {}'.format(csw.sw)) #don't use tretrv, it doesn't work

    iyd,utsec,stl = datetime2gtd(dtime,glon)

    altkm = atleast_1d(altkm)
    dens = empty((altkm.size,9)); temp=empty((altkm.size,2))

    meters(1) # makes output in m^-3 and kg/m^-3
    for i,a in enumerate(altkm):
        dens[i,:],temp[i,:] = gtd7(iyd,utsec,a,glat,glon,stl, f107a,f107, ap,mass)

    densd = DataFrame(dens, index=altkm, columns=species)
    tempd = DataFrame(temp, index=altkm, columns=ttypes)
    return densd,tempd


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

    dens,temp = rungtd1d(dtime,altkm,glat,glon,p.f107a,p.f107,p.ap,p.mass)
