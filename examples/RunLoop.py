#!/usr/bin/env python3
"""
To generate eigenprofiles, use -E option.

Poker Flat eigenprofile example from command line:
python 2013-03-01T10:48Z -c 65 -148 -E ~/code/transcar/transcar/BT_E1E2prev.csv --f107 --f107a --ap

Example mar 1 2011:
------------------
1) create unit input flux spectrum with gridaurora/MakeEigenprofileFluxInput.py
2) run this program:
python3 RunLoop.py -t 2011-03-01T00:00Z -E ~/data/100MeVtop.h5 -c 65 -147.5 --f107 115 --f107p 115 --f107a 96 --ap 7 -o ~/data/rates.h5



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
from collections import namedtuple
from datetime import datetime
from pytz import UTC
from dateutil import rrule
from dateutil.parser import parse
from matplotlib.pyplot import show
#
from glowaurora.runglow import plotprodloss,plotaurora
from glowaurora.eigenprof import verprodloss,makeeigen,ekpcolor
from histfeas.plotsnew import ploteigver

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="Stan Solomon's GLOW auroral model")
    p.add_argument('-t','--simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',nargs='+',required=True)#,default='1999-12-21T00:00:00Z')
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(70,0))
    p.add_argument('--flux',help='overall incident flux [erg ...]',type=float,default=1.)
    #grp = parser.add_mutually_exclusive_group()

    #'~/code/transcar/transcar/BT_E1E2prev.csv')  #~/data/100MeVtop.h5
    p.add_argument('-E','--eigenprof',help='generate eigenprofiles using energies in this HDF5 or CSV file')
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
        T = [parse(p.simtime[0])]
    elif len(p.simtime) == 2:
        T = list(rrule.rrule(rrule.HOURLY,
                                 dtstart=parse(p.simtime[0]),
                                 until =parse(p.simtime[1])))

    makeplot = p.makeplot

    if p.eigenprof: #loop over time
        makeplot.append('eig')
        EKpcolor,EK,diffnumflux = ekpcolor(p.eigenprof)

        ver,photIon,isr,phitop,zceta,sza,EKpcolor,prates,lrates,sion = makeeigen(EK,diffnumflux,T,p.latlon,
                                                                             p.f107a,p.f107,p.f107p,p.ap,
                                                                             p.makeplot,p.odir,p.zlim)
    else: #single time
        ver,photIon,isr,phitop,zceta,sza,prates,lrates,sion = verprodloss(T[0],p.latlon,p.flux,p.e0,
                                                                       p.f107a,p.f107,p.f107p,p.ap,
                                                                       p.makeplot,p.odir,p.zlim)

#%% plotting
    if p.eigenprof:
        z=ver.major_axis.values
        sim = namedtuple('sim',['reacreq','opticalfilter']); sim.reacreq=sim.opticalfilter=''
        zlim=(None,None)
        glat=p.latlon[0]; glon=p.latlon[1]
        z=ver.major_axis.values
        for t in ver: #for each time
            #VER eigenprofiles
            ploteigver(EKpcolor,z,ver[t].values,(None,)*6,sim,str(t)+' Vol. Emis. Rate ')

            if False:
                for prate,lrate,E0 in zip(prates,lrates,EKpcolor):
                    #production eigenprofiles
                    plotprodloss(z,prate,lrate,t,glat,glon,zlim,'Volume Production/Loss Rates',' E0: {:.0f}'.format(E0))
                    #loss eigenprofiles
    else:
        plotaurora(phitop,ver,zceta,photIon,isr,sion,T[0],p.latlon[0],p.latlon[1],prates,lrates,makeplot=makeplot)

    show()
