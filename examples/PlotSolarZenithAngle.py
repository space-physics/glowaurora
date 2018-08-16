#!/usr/bin/env python3
"""
Comparing solzen.f solar zenith angle with Astropy, per email discussion with UAF and UCAR staff
Michael Hirsch
2015
"""
import numpy as np
from pandas import date_range
import xarray
from matplotlib.pyplot import show, figure
from time import time
from datetime import datetime
from dateutil.parser import parse
import astropy.units as u
from astropy.coordinates import get_sun, EarthLocation, AltAz
from astropy.time import Time
#
from sciencedates import datetime2yd
#
import glow as glowfort


def demosolzen(t: datetime, glat: float, glon: float):
    """demo the solar zenith angle calclation vs AstroPy"""
# %% SZA with glow
    yd, utsec = datetime2yd(dtime)

    sza_glow = np.empty_like(dtime, dtype=float)
    for j, (d, s) in enumerate(zip(yd, utsec)):
        sza_glow[j] = glowfort.sun_angles.solzen(d, s, glat, glon)

    return sza_glow


def demosuncor(T, glat, glon, alt_m):
    # %% Solar location with GLOW
    yd, utsec = datetime2yd(T)[:2]
    solar = xarray.DataArray(np.empty((T.size, 2, 4)),
                             dims=['time', 'algorithm', 'coord'],
                             coords={'time': T,
                                     'algorithm': ['glow', 'astropy'],
                                     'coord': ['dec', 'ra', 'gst', 'sza']})
    tic = time()
    for d, s, t in zip(yd, utsec, T):
        solar.loc[t, 'glow', ['dec', 'ra', 'gst']] = glowfort.sun_angles.suncor(d, s)
    fortsec = time() - tic
    solar.loc[:, 'glow', ['dec', 'ra', 'gst']] = np.degrees(solar.loc[:, 'glow', ['dec', 'ra', 'gst']])
# %% solar location with AstroPy
    tic = time()
    sza, sun, sunobs = solarzenithangle(T, glat, glon, alt_m)
    pysec = time() - tic
    solar.loc[:, 'astropy', 'dec'] = sun.dec.degree
    solar.loc[:, 'astropy', 'ra'] = sun.ra.degree
    solar.loc[:, 'astropy', 'gst'] = sunobs.obstime.sidereal_time('apparent', 'greenwich').degree
    solar.loc[:, 'astropy', 'sza'] = sza

    print(f'in seconds, fortran time: {fortsec:.3f}   python time: {pysec:.3f} ')

    return solar


def plotsza(sun, error):
    ax = figure().gca()
    ax.plot(sun.time, sun.loc[:, 'glow', 'sza'], label='Glow')
    ax.plot(sun.time, sun.loc[:, 'astropy', 'sza'], label='AstroPy')
    ax.set_title('comparison of solar zenith angle [degrees]')

    ax = figure().gca()
    ax.plot(error.time, error.loc[:, 'sza'])
    ax.set_title('SINGLE PRECISION Error in GLOW solar predication [degrees]')


def solarzenithangle(t, glat, glon, alt_m):
    """
    Input:

    t: scalar or array of datetime

    """
    if isinstance(t, str):
        t = [t]
    if not isinstance(t[0], datetime):
        t = map(parse, t)

    obs = EarthLocation(lat=glat * u.deg, lon=glon * u.deg, height=alt_m * u.m)
    times = Time(t, scale='ut1')
    sun = get_sun(times)
    sunobs = sun.transform_to(AltAz(obstime=times, location=obs))

    return 90. - sunobs.alt.degree, sun, sunobs


if __name__ == '__main__':
    dtime = date_range('2013-04-01', '2013-04-02',
                       freq='2Min', tz='utc', closed='left').to_pydatetime()
    lat = 65
    lon = -148
    alt_m = 0  # to be fair to GLOW
# %% compute coordinates of sun
    sun = demosuncor(dtime, lat, lon, alt_m)
# %% compute angle from zenith of sun
    sun.loc[:, 'glow', 'sza'] = demosolzen(dtime, lat, lon)
# %% compute error
    error = sun.loc[:, 'astropy', :] - sun.loc[:, 'glow', :]

    plotsza(sun, error)
    show()

    assert error.max() < 0.005
