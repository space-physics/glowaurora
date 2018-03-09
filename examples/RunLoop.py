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
from dateutil import rrule
from dateutil.parser import parse
from matplotlib.pyplot import show
#
from glowaurora.plots import plotprodloss,plotaurora
from glowaurora import verprodloss,makeeigen,ekpcolor
try:
    from histfeas.plotsnew import ploteigver
except ImportError:
    ploteigvar=None

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="Stan Solomon's GLOW auroral model")
    p.add_argument('-t','--simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',nargs='+',default=['2013-04-14T15:54Z'])
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(65,-148))
    p.add_argument('--flux',help='overall incident flux [erg ...]',type=float,default=1.)
    #grp = parser.add_mutually_exclusive_group()

    #'~/code/transcar/transcar/BT_E1E2prev.csv')  #~/data/100MeVtop.h5
    p.add_argument('-E','--eigenprof',help='generate eigenprofiles using energies in this HDF5 or CSV file')
    p.add_argument('--e0',help='characteristic energy [eV]',type=float,nargs='+',default=(1e3,))

    p.add_argument('-m','--makeplot',help='show to show plots, png to save pngs of plots',nargs='+',default=['show'])
    p.add_argument('-o','--odir',help='output directory to write plots',default='')
    p.add_argument('-z','--zlim',help='minimum,maximum altitude [km] to plot',nargs=2,default=(60,350),type=float)
    p = p.parse_args()

    params = {'glat':p.latlon[0],
              'glon':p.latlon[1],
              'flux':p.flux,
              'EK':p.e0,
              'makeplot':p.makeplot,
              'zlim':p.zlim,
              'plotformat':'png',
            }

    if len(p.simtime) == 1:
        params['t0'] = [parse(p.simtime[0])]
    elif len(p.simtime) == 2:
        params['t0'] = list(rrule.rrule(rrule.HOURLY,
                                 dtstart=parse(p.simtime[0]),
                                 until =parse(p.simtime[1])))

    if p.eigenprof: #loop over time
        params['makeplot'].append('eig')
        EKpcolor,EK,diffnumflux = ekpcolor(p.eigenprof)

        sim = makeeigen(params)
    else: #single time
        sim = verprodloss(params)

#%% plotting
    if p.eigenprof:
        z=sim['ver'].major_axis.values
        optsim = namedtuple('sim',['reacreq','opticalfilter']); sim.reacreq=sim.opticalfilter=''
        zlim=(None,None)
        glat=p.latlon[0]; glon=p.latlon[1]
        z = sim.z_km
        for t in sim.time: #for each time
            #VER eigenprofiles
            if ploteigver:
                ploteigver(EKpcolor,z,sim['ver'][t].values, (None,)*6, optsim,
                           f'{t} Vol. Emis. Rate ')
            else:
                for prate,lrate,E0 in zip(sim['prates'],sim['lrates'],EKpcolor):
                    #production eigenprofiles
                    plotprodloss(z,prate,lrate,t,glat,glon,zlim,
                                 f'Volume Production/Loss Rates',' E0: {E0:.0f}')
                    #loss eigenprofiles
    else:
        plotaurora(params, sim)

    show()
