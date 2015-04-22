#!/usr/bin/env python3
"""
Trivial example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
bostonmicrowave.com
"""
from matplotlib.pyplot import figure, show,subplots
from pandas import DataFrame
from datetime import datetime
from dateutil.parser import parse
from numpy import hstack
#
from fortrandates import datetime2gtd
try:
    from aurora import aurora
except ImportError as e:
    exit('you must compile with f2py first. See README.md  {}'.format(e))

def demoaurora(iyd,utsec,glat,glon,f107a,f107,f107p,ap):
    z,zeta,ion,ecalc,photI,ImpI = aurora(iyd,utsec,glat,glon,f107a,f107,f107p,ap,1,1000)

    ver = DataFrame(index=z,
                    data=zeta.T[:,:10],
                    columns=[3371, 4278, 5200, 5577, 6300,7320,10400,3466,
                             7774, 8446])
    photIon = DataFrame(index=z,
                   data=hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)),
                    columns=['photoIoniz','eImpactIoniz','eDens',
                    'nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO'])
    return ver,photIon

def plotaurora(ver,photIon,dtime,glat,glon):
    fg,axs = subplots(1,3,sharey=True)

    ax = axs[0]
    ax.plot(ver.values,ver.index)
    ax.set_xlabel('VER')
    ax.set_ylabel('altitude [km]')
    ax.grid(True)
    ax.legend(ver.columns)
    ax.set_title('{}  ({},{})'.format(dtime,glat,glon))

    ax = axs[1]
    ax.plot(photIon[['photoIoniz','eImpactIoniz']],photIon.index)
    ax.set_xlabel('ionization')
    ax.grid(True)
    ax.legend(photIon.columns)
    ax.set_title('{}  ({},{})'.format(dtime,glat,glon))

    ax = axs[2]
    ax.semilogx(photIon[['nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO']], photIon.index)
    ax.set_xlabel('density')
    ax.grid(True)
    ax.legend(photIon.columns)
    ax.set_title('{}  ({},{})'.format(dtime,glat,glon))

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="Stan Solomon's GLOW auroral model")
    p.add_argument('simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',type=str,nargs='?',default=None)
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(65,-148))
    p.add_argument('--f107a',help='AVERAGE OF F10.7 FLUX',type=float,default=150)
    p.add_argument('--f107p',help='DAILY F10.7 FLUX FOR PREVIOUS DAY',type=float,default=150)
    p.add_argument('--f107',help='F10.7 for sim. day',type=float,default=150)
    p.add_argument('--ap',help='daily ap',type=float,default=4)
    p = p.parse_args()

    if p.simtime is None:
        dtime = datetime.now()
    else:
        dtime = parse(p.simtime)
    iyd,utsec = datetime2gtd(dtime)[:2]

    (glat,glon) = p.latlon

    ver,photIon = demoaurora(iyd,utsec,glat,glon,p.f107a,p.f107,p.f107p,p.ap)
    plotaurora(ver,photIon,dtime,glat,glon)
    show()
