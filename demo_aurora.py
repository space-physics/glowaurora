#!/usr/bin/env python3
"""
Trivial example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
bostonmicrowave.com
"""
from matplotlib.pyplot import figure, show,subplots,tight_layout
from pandas import DataFrame
from datetime import datetime
from dateutil.parser import parse
from numpy import hstack,arange,append,array
try:
    import seaborn as sns
except ImportError as e:
    print('Seaborn not installed, falling back to basic Matplotlib plots.  {}'.format(e))
#
from fortrandates import datetime2gtd
try:
    from glowgrid import energygrid,maxt
    from aurora import aurora
except ImportError as e:
    exit('you must compile with f2py first. See README.md  {}'.format(e))

def demoaurora(nbins,eflux,e0,iyd,utsec,glat,glon,f107a,f107,f107p,ap):
#%% temporarily use glow grid instead of our own
    ener,dE = energygrid(nbins)
    phitop = maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)
    phi = hstack((ener[:,None],dE[:,None],phitop[:,None]))
#%% glow model
    z = arange(80,110+1,1)
    z = append(z,array([111.5,113.,114.5,116.,118.,120.,122.,124.,126., 128.,130.,132.,134.,136.,138.,140.,142.,144.,146., 148.,150.,153.,156.,159.,162.,165.,168.,172.,176., 180.,185.,190.,195.,200.,205.,211.,217.,223.,230.,237.,244.,252.,260.,268.,276.,284.,292.,300.,309., 318.,327.,336.,345.,355.,365.,375.,385.,395.,406., 417.,428.,440.,453.,467.,482.,498.,515.,533.,551., 570.,590.,610.,630.,650.,670.,690.,710.,730.,750., 770.,790.,810.,830.,850.,870.,890.,910.,930.,950.]))

    zeta,ion,ecalc,photI,ImpI,isr = aurora(z,iyd,utsec,glat,glon%360,
                                             f107a,f107,f107p,ap,phi)

    ver = DataFrame(index=z,
                    data=zeta[:,:10],
                    columns=[3371, 4278, 5200, 5577, 6300,7320,10400,3466,
                             7774, 8446])
    photIon = DataFrame(index=z,
                   data=hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)),
                    columns=['photoIoniz','eImpactIoniz','ne',
                    'nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO'])

    isrparam = DataFrame(index=z,
                         data=isr,
                         columns=['ne','Te','Ti'])

    phitop = DataFrame(index=phi[:,0],
                       data=phi[:,2],
                       columns=['flux'])
    return ver,photIon,isrparam,phitop

def plotaurora(phitop,ver,photIon,isr,dtime,glat,glon):
    ax = figure().gca()
    phitop.plot(ax=ax,logx=True,logy=True)
    ax.set_title('Incident Flux',fontsize='x-large')
    ax.set_xlabel('Beam Energy [eV]',fontsize='large')
    ax.set_ylabel('Flux',fontsize='large')
    ax.tick_params(axis='both',which='major',labelsize='medium')


    fg,axs = subplots(1,4,sharey=True, figsize=(15,8))
    fg.suptitle('{} ({},{})'.format(dtime,glat,glon),fontsize='x-large')
    tight_layout(pad=3.2, w_pad=0.3)

    ax = axs[0]
    ax.plot(ver.values,ver.index)
    ax.set_xlabel('VER',fontsize='large')
    ax.set_ylabel('altitude [km]',fontsize='large')
    ax.legend(ver.columns)
    ax.set_title('Volume emission rate',fontsize='x-large')

    ax = axs[1]
    ax.plot(photIon[['photoIoniz','eImpactIoniz']],photIon.index)
    ax.set_xlabel('ionization',fontsize='large')
    ax.legend(photIon.columns[:2])
    ax.set_title('Photo and e$^-$ impact ionization',fontsize='x-large')

    ax = axs[2]
    ax.semilogx(photIon[['ne','nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO']], photIon.index)
    ax.set_xlabel('Density',fontsize='large')
    ax.legend(photIon.columns[2:])
    ax.set_title('Electron and Ion Densities',fontsize='x-large')

    ax = axs[3]
    ax.semilogx(isr[['Te','Ti']], isr.index)
    ax.set_xlabel('Temperature [K]',fontsize='large')
    ax.legend(isr.columns[1:])
    ax.set_title('Particle Temperature',fontsize='x-large')

    for a in axs:
        a.grid(True)
        a.tick_params(axis='both',which='major',labelsize='medium')

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
    iyd,utsec = datetime2gtd(dtime)[:2]

    (glat,glon) = p.latlon

    ver,photIon,isr,phitop = demoaurora(p.nbins,p.flux,p.e0,iyd,utsec,glat,glon,
                                        p.f107a,p.f107,p.f107p,p.ap)
    plotaurora(phitop,ver,photIon,isr,dtime,glat,glon)
    show()
