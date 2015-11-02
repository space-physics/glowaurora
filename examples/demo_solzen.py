#!/usr/bin/env python3
"""
Comparing solzen.f solar zenith angle with Astropy, per email discussion with UAF and UCAR staff
Michael Hirsch
2015
"""
from __future__ import division,absolute_import
from six.moves import getcwd
from numpy import empty_like,degrees
from pandas import date_range,DataFrame,Panel
from matplotlib.pyplot import show
from os import chdir,environ
from time import time
try:
    import seaborn
except:
    pass
#
from histutils.fortrandates import datetime2yd
from gridaurora.solarangle import solarzenithangle
#################################
#TODO hack for module data path issue
chdir(environ['HOME'])
import glowaurora
from glowaurora import glowfort
chdir(glowaurora.__path__[0])
print('loaded glow from ' + getcwd())
#################################
#%% demo the solar zenith angle calclation vs AstroPy
def demosolzen(dtime,glat,glon):
#%% SZA with glow
    yd,utsec = datetime2yd(dtime)[:2]

    sza_glow = empty_like(dtime,dtype=float)
    for j,(d,s) in enumerate(zip(yd,utsec)):
        sza_glow[j] = glowfort.solzen(d,s,glat,glon)

    return DataFrame(index=dtime,data=sza_glow,columns=['glow'])

def demosuncor(T,glat,glon,alt_m):
#%% Solar location with GLOW
    yd,utsec = datetime2yd(T)[:2]
    solar = Panel( items=['glow','astropy'],
                 major_axis=T,
                 minor_axis=['dec','ra','gst'])
    tic = time()
    for d,s,t in zip(yd,utsec,T):
        solar['glow'].loc[t] = glowfort.suncor(d,s)
    fortsec = time()-tic
    solar.ix['glow',:,'dec'] = degrees(solar.ix['glow',:,'dec'])
    solar.ix['glow',:,'ra'] = degrees(solar.ix['glow',:,'ra'])
    solar.ix['glow',:,'gst'] = degrees(solar.ix['glow',:,'gst'])
#%% solar location with AstroPy
    tic=time()
    sza,sun,sunobs = solarzenithangle(T,glat,glon,alt_m)
    pysec = time()-tic
    solar.ix['astropy',:,'dec'] = sun.dec.degree
    solar.ix['astropy',:,'ra'] = sun.ra.degree
    solar.ix['astropy',:,'gst'] = sunobs.obstime.sidereal_time('apparent','greenwich').degree

    print('in seconds, fortran time: {:.3f}   python time: {:.3f} '.format(fortsec,pysec))

    return solar,sza



def plotsza(sza,sun,error):
    sza.plot(title='comparison of solar zenith angle [degrees]')

    error.plot(subplots=True,title='SINGLE PRECISION Error in GLOW solar predication [degrees]')

if __name__ == '__main__':
    dtime = date_range('2013-04-01','2013-04-02',
                       freq='2Min',tz='utc',closed='left').to_pydatetime()
    lat=65
    lon=-148
    alt_m = 0 #to be fair to GLOW
#%% compute coordinates of sun
    sun,szaastropy = demosuncor(dtime,lat,lon,alt_m)
#%% compute angle from zenith of sun
    sza = demosolzen(dtime,lat,lon)
    sza['astropy'] = szaastropy
#%% compute error
    error = sun['astropy'] - sun['glow']
    error['sza'] = sza['astropy'] - sza['glow']

    plotsza(sza,sun,error)
    show()

    assert (error.max() < [0.35,0.75,0.005,0.003]).all()
    print('OK')