from pathlib import Path
import logging
from numpy import repeat,empty
from xarray import DataArray
from numpy import atleast_1d
from os import chdir
#
from sciencedates import datetime2gtd
#
import glowaurora
from glowfort import tselec,meters,csw,gtd7
glowpath=Path(glowaurora.__path__[0])/'..'
oldcwd = Path.cwd()

def rungtdGLOW(dtime,altkm,glat,glon,f107a,f107,ap,mass,tselecopts):
    chdir(str(glowpath))
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

    densd = DataArray(dens, dims=['z_km','species'],
                      coords={'z_km':altkm, 'species':species})
    tempd = DataArray(temp, dims=['z_km','temperature'],
                      coords={'z_km':altkm, 'temperature':ttypes})
    return densd,tempd
