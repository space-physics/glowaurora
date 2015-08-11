#!/usr/bin/env python3
from datetime import datetime
from dateutil.parser import parse
from matplotlib.pyplot import show
#
from gridaurora.fortrandates import datetime2yd
from glowaurora.runglow import glowaurora,plotaurora

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="Stan Solomon's GLOW auroral model")
    p.add_argument('simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',type=str,nargs='?',default=None)
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(65,-148))
    p.add_argument('-n','--nbins',help='number of energy bins in incident diff num flux',type=int,default=190)
    p.add_argument('--flux',help='overall incident flux [erg ...]',type=float,default=1)
    p.add_argument('--e0',help='characteristic energy [eV]',type=float,default=1e3)
    p.add_argument('--f107a',help='AVERAGE OF F10.7 FLUX',type=float,default=150)
    p.add_argument('--f107p',help='DAILY F10.7 FLUX FOR PREVIOUS DAY',type=float,default=150)
    p.add_argument('--f107',help='F10.7 for sim. day',type=float,default=150)
    p.add_argument('--ap',help='daily ap',type=float,default=4)
    p = p.parse_args()

    if p.simtime is None:
        dtime = datetime.now()
    else:
        dtime = parse(p.simtime)

    yd,utsec = datetime2yd(dtime)[:2]

    (glat,glon) = p.latlon

    ver,photIon,isr,phitop,zceta = glowaurora(p.nbins,p.flux,p.e0,
                                              yd,utsec,glat,glon,
                                        p.f107a,p.f107,p.f107p,p.ap)

    plotaurora(phitop,ver,zceta,photIon,isr,dtime,glat,glon)
    show()