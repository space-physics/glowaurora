from numpy import empty_like,degrees
from pandas import date_range,DataFrame,Panel
from matplotlib.pyplot import figure,show
# Some features may require AstroPy 1.0+
import astropy.units as u
from astropy.coordinates import get_sun, EarthLocation, AltAz
from astropy.time import Time
#
from gridaurora.fortrandates import datetime2yd
from aurora import szacalc
#%% demo the solar zenith angle calclation vs AstroPy
def demosolzen(dtime,glat,glon):
#%% SZA with glow
    yd,utsec = datetime2yd(dtime)[:2]

    sza_glow = empty_like(dtime,dtype=float)
    for j,(d,s) in enumerate(zip(yd,utsec)):
        sza_glow[j] = szacalc.solzen(d,s,glat,glon)

    return DataFrame(index=dtime,data=sza_glow,columns=['glow'])

def demosuncor(dtime,glat,glon,alt_m):
#%% Solar location with GLOW
    yd,utsec = datetime2yd(dtime)[:2]
    solar = Panel( items=['glow','astropy'],
                 major_axis=dtime,
                 minor_axis=['dec','ra','gst'])
    for d,s,t in zip(yd,utsec,dtime):
        solar['glow'].loc[t] = szacalc.suncor(d,s)
    solar.ix['glow',:,'dec'] = degrees(solar.ix['glow',:,'dec'])
    solar.ix['glow',:,'ra'] = degrees(solar.ix['glow',:,'ra'])
    solar.ix['glow',:,'gst'] = degrees(solar.ix['glow',:,'gst'])
#%% solar location with AstroPy
    obs = EarthLocation(lat=glat*u.deg, lon=glon*u.deg, height=alt_m*u.m)
    times = Time(dtime, scale='ut1')
    sun = get_sun(times)
    sunobs = sun.transform_to(AltAz(obstime=times,location=obs))
    solar.ix['astropy',:,'dec'] = sun.dec.degree
    solar.ix['astropy',:,'ra'] = sun.ra.degree
    solar.ix['astropy',:,'gst'] = sunobs.obstime.sidereal_time('apparent','greenwich').degree

    return solar,sunobs



def plotsza(sza,sun,error):
    ax=figure().gca()
    sza.plot(ax=ax,subplots=True,title='comparison of solar zenith angle [degrees]')

    ax = figure().gca()
    error.plot(ax=ax,subplots=True,title='1973 SINGLE PRECISION Error in GLOW solar predication [degrees]')

if __name__ == '__main__':
    dtime = date_range('2013-04-01','2013-04-02',
                       freq='2Min',tz='utc',closed='left').to_pydatetime()
    lat=65
    lon=-148
    alt_m = 0 #to be fair to GLOW
#%% compute coordinates of sun
    sun,sunobs = demosuncor(dtime,lat,lon,alt_m)
#%% compute angle from zenith of sun
    sza = demosolzen(dtime,lat,lon)
    sza['astropy'] = 90 - sunobs.alt.degree
#%%
    error = sun['astropy'] - sun['glow']
    error['sza'] = sza['astropy'] - sza['glow']
    plotsza(sza,sun,error)
    show()
