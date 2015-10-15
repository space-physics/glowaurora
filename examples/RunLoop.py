#!/usr/bin/env python3
"""
To generate eigenprofiles, use -E option.

Poker Flat eigenprofile example from command line:
python 2013-03-01T10:48Z -c 65 -148 -E ~/code/transcar/transcar/BT_E1E2prev.csv --f107 --f107a --ap

example mar 1 2011
python3 RunLoop.py -t 2011-03-01T00:00Z 2013-03-01T23:00Z -E ~/code/transcar/transcar/BT_E1E2prev.csv -c 65 -148 --f107 115 --f107p 115 --f107a 96 --ap 7

default parameter values like those of Stan's fortran examples--yield rather similar output
Note that the number of bins for altitude and energy are "compiled in" the Fortran.
It would require substantial changes to the Fortran code to make these dynamic due to Common Block
usage in the Fortran77 code. Stan Solomon is upgrading the GLOW code to Fortran 90 with modules,
so I will perhaps try again with dynamic bin size for altitude and energy when
that code is available.

f10.7 and ap available from ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt

Michael Hirsch
Sept 2015
"""
from __future__ import division,absolute_import
import h5py
from datetime import datetime
from pytz import UTC
from dateutil import rrule
from dateutil.parser import parse
from matplotlib.pyplot import show
from os.path import expanduser
from numpy import loadtxt,append,asarray
from pandas import DataFrame,Panel
#
from glowaurora.runglow import runglowaurora#,plotaurora
from histfeas.plotsnew import ploteig
from transcarread.readTranscar import SimpleSim

epoch = datetime(1970,1,1,tzinfo=UTC)

def E0aurora(dt,glatlon,flux,E0,f107a,f107,f107p,ap,makeplot,odir,zlim):

    (glat,glon) = glatlon

    DFver = DataFrame()
    for e0 in E0:
        print('char. energy {}'.format(e0))

        ver,photIon,isr,phitop,zceta,sza = runglowaurora(flux,e0,
                                              dt,glat,glon,
                                              f107a,f107,f107p,ap)

        #plotaurora(phitop,ver,flux,sza,zceta,photIon,isr,dtime,glat,glon,e0,zlim,makeplot,odir)

        DFver[e0] = ver.sum(axis=1)

    return DFver,photIon,isr,phitop,zceta,sza

def ekpcolor(eigenfn):
    e0 =   loadtxt(expanduser(eigenfn),usecols=[0],delimiter=',')
    eEnd = loadtxt(expanduser(eigenfn),usecols=[1],delimiter=',')[-1]
    return append(e0,eEnd),e0

def makeeigen(eigenfn,dt,glatlon,f107a,f107,f107p,ap,makeplot,odir,zlim):
    makeplot.append('eig')
    flux=None
    EKpcolor,e0 = ekpcolor(eigenfn)

    ver = None

    for t in dt:
        v,photIon,isr,phitop,zceta,sza = E0aurora(dtime,p.latlon,flux,e0,
                                            p.f107a,p.f107,p.f107p,p.ap,
                                            p.makeplot,p.odir,p.zlim)
        if ver is None:
            ver = Panel(items=dt,major_axis=v.index,minor_axis=v.columns)
        ver.loc[t,:,:] = v

    return ver,photIon,isr,phitop,zceta,sza,EKpcolor

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="Stan Solomon's GLOW auroral model")
    p.add_argument('-t','--simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',nargs='+',required=True)#,default='1999-12-21T00:00:00Z')
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(70,0))
    p.add_argument('--flux',help='overall incident flux [erg ...]',type=float,default=1.)
    #grp = parser.add_mutually_exclusive_group()
    p.add_argument('-E','--eigenprof',help='generate eigenprofiles using energies in this csv file')#'~/code/transcar/transcar/BT_E1E2prev.csv')
    p.add_argument('--e0',help='characteristic energy [eV]',type=float,nargs='+',default=(1e3,))

    p.add_argument('--f107a',help='AVERAGE OF F10.7 FLUX',type=float,default=100)
    p.add_argument('--f107p',help='DAILY F10.7 FLUX FOR PREVIOUS DAY',type=float,default=100)
    p.add_argument('--f107',help='F10.7 for sim. day',type=float,default=100)
    p.add_argument('--ap',help='daily ap',type=float,default=4)
    p.add_argument('-m','--makeplot',help='show to show plots, png to save pngs of plots',nargs='+',default=['show'])
    p.add_argument('-o','--odir',help='output directory to write plots',default='')
    p.add_argument('-z','--zlim',help='minimum,maximum altitude [km] to plot',nargs=2,default=(60,350),type=float)
    p = p.parse_args()

    if len(p.simtime) == 1:
        dtime = [parse(p.simtime[0])]
    elif len(p.simtime) == 2:
        dtime = list(rrule.rrule(rrule.HOURLY,
                                 dtstart=parse(p.simtime[0]),
                                 until =parse(p.simtime[1])))

    makeplot = p.makeplot

    if p.eigenprof:
        ver,photIon,isr,phitop,zceta,sza,EKpcolor = makeeigen(p.eigenprof,dtime,p.latlon,
                                                             p.f107a,p.f107,p.f107p,p.ap,
                                                             p.makeplot,p.odir,p.zlim)
    else:
        flux=p.flux; e0=p.e0

        ver,photIon,isr,phitop,zceta,sza = E0aurora(dtime,p.latlon,flux,e0,
                                                p.f107a,p.f107,p.f107p,p.ap,
                                                p.makeplot,p.odir,p.zlim)

    if p.odir.endswith('.h5'):
        h5fn = expanduser(p.odir)
        print('writing to '+h5fn)
        ut1_unix = [(t-epoch).total_seconds() for t in ver.items.to_pydatetime()]
        with h5py.File(h5fn,'w',libver='latest') as f:
            d=f.create_dataset('/eigenprofile',data=ver.values)
            d=f.create_dataset('/altitude',data=ver.major_axis)
            d=f.create_dataset('/Ebins',data=ver.minor_axis)
            d=f.create_dataset('/ut1_unix',data=ut1_unix)
#%% plotting
    if p.eigenprof:
        sim = SimpleSim(filt='none',inpath=None,reacreq='')
        for t in ver: #for each time
            ploteig(EKpcolor,asarray(ver.major_axis),ver[t],(None,)*6,sim,t)
        show()

