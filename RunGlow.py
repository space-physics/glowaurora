#!/usr/bin/env python
"""
default parameter values like those of Stan's fortran examples--yield rather similar output

Note, this simulation uses a specific input differential number flux spectrum
"""
from matplotlib.pyplot import show
#import seaborn
#
import glowaurora as glow
from glowaurora.plots import plotaurora

def E0aurora(t0,glatlon,flux,e0,f107a,f107,f107p,ap,makeplot):

    (glat,glon) = glatlon

    ver,photIon,isr,phitop,zceta,sza,prate,lrate,tez,sion = glow.runglowaurora(flux,e0,
                                                                 t0,glat,glon,
                                                                 f107a,f107,f107p,ap)

    plotaurora(phitop,ver,zceta,photIon,isr,sion,t0,glat,glon,prate,lrate,tez,e0,makeplot=makeplot)

    return ver,photIon,isr,phitop,zceta,sza

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="Stan Solomon's GLOW auroral model")
    p.add_argument('-t','--simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',default='2013-04-14T15:54Z')
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(65.,-148.))
#    p.add_argument('-n','--nbins',help='number of energy bins in incident diff num flux',type=int,default=190) #hard-coded in cglow.h
    p.add_argument('-q','--flux',help='overall incident flux [erg ...]',type=float,default=1.)
    p.add_argument('--e0',help='characteristic energy [eV]',type=float,default=1e3)
    p.add_argument('--f107a',help='AVERAGE OF F10.7 FLUX',type=float)
    p.add_argument('--f107p',help='DAILY F10.7 FLUX FOR PREVIOUS DAY',type=float)
    p.add_argument('--f107',help='F10.7 for sim. day',type=float)
    p.add_argument('--ap',help='daily ap',type=float)
    p.add_argument('-m','--makeplot',help='show to show plots, png to save pngs of plots',nargs='+',default=['show'])
    p = p.parse_args()

    ver,photIon,isr,phitop,zceta,sza = E0aurora(p.simtime,p.latlon,p.flux,p.e0,p.f107a,p.f107,p.f107p,p.ap,p.makeplot)

    show()
