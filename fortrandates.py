from six import string_types
from datetime import timedelta,datetime, time
from pytz import UTC
from numpy import atleast_1d, empty_like, atleast_2d,nan,empty
from dateutil.parser import parse

def datetime2yd(dtime):
    """
    Inputs:
    dtime: Numpy 1-D array of datetime.datetime OR string suitable for dateutil.parser.parse

    Outputs:
    yd: yyyyddd four digit year, 3 digit day of year (INTEGER)
    utsec: seconds from midnight utc
    """
    dtime = atleast_1d(dtime)

    utsec=empty_like(dtime,dtype=float)
    yd = empty_like(dtime,dtype=int)
    for i,t in enumerate(dtime):
        if isinstance(t,string_types):
            t = parse(t)

        t=forceutc(t)
        utsec[i] = dt2utsec(t)
        yd[i] = t.year*1000 + int(t.strftime('%j'))

    return yd,utsec


def datetime2gtd(dtime,glon=nan):
    """
    Inputs:
    dtime: Numpy 1-D array of datetime.datetime OR string suitable for dateutil.parser.parse
    glon: Numpy 2-D array of geodetic longitudes

    Outputs:
    doy: day of year
    utsec: seconds from midnight utc
    stl: local solar time
    """
    dtime = atleast_1d(dtime); glon=atleast_2d(glon)
    doy=empty_like(dtime,dtype=int); utsec=empty_like(dtime,dtype=float)
    stl = empty((dtime.size,glon.shape[0],glon.shape[1]))

    for i,t in enumerate(dtime):
        if isinstance(t,string_types):
            t = parse(t)

        t = forceutc(t)

        doy[i] = int(t.strftime('%j'))
        utsec[i] = dt2utsec(t)

        stl[i,...] = utsec[i]/3600 + glon/15 #FIXME let's be sure this is appropriate

    return doy,utsec,stl

def dt2utsec(dt):
    """ seconds since utc midnight"""
    return timedelta.total_seconds(dt-datetime.combine(dt.date(),time(0,tzinfo=UTC)))


def forceutc(t):
    """
    input: python datetime (naive, utc, non-utc)
    output: utc datetime
    """
    if t.tzinfo == None:
        t = t.replace(tzinfo = UTC)
    else:
        t = t.astimezone(UTC)
    return t

"""
http://stackoverflow.com/questions/19305991/convert-fractional-years-to-a-real-date-in-python
Authored by "unutbu" http://stackoverflow.com/users/190597/unutbu

In Python, go from decimal year (YYYY.YYY) to datetime,
and from datetime to decimal year.
"""
def yeardec2datetime(atime):
    """
    Convert atime (a float) to DT.datetime
    This is the inverse of dt2t.
    assert dt2t(t2dt(atime)) == atime
    """
    year = int(atime)
    remainder = atime - year
    boy = datetime(year, 1, 1)
    eoy = datetime(year + 1, 1, 1)
    seconds = remainder * (eoy - boy).total_seconds()
    return boy + datetime.timedelta(seconds=seconds)

def datetime2yeardec(adatetime):
    """
    Convert a datetime into a float. The integer part of the float should
    represent the year.
    Order should be preserved. If adate<bdate, then d2t(adate)<d2t(bdate)
    time distances should be preserved: If bdate-adate=ddate-cdate then
    dt2t(bdate)-dt2t(adate) = dt2t(ddate)-dt2t(cdate)
    """
    if isinstance(adatetime,string_types):
        adatetime = parse(adatetime)

    year = adatetime.year
    boy = datetime(year, 1, 1)
    eoy = datetime(year + 1, 1, 1)
    return year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds()))
