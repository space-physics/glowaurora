#!/usr/bin/env python3
"""
Trivial example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
"""
from __future__ import division,absolute_import
from itertools import chain
from matplotlib.pyplot import figure, subplots,tight_layout
from pandas import DataFrame
from numpy import hstack,asarray,rollaxis
from os import chdir
try:
    import seaborn
except:
    pass
#
from histutils.fortrandates import datetime2yd
import glowaurora
from glowaurora import glowfort
#
glowpath=glowaurora.__path__[0]

def runglowaurora(eflux,e0,dt,glat,glon,f107a,f107,f107p,ap):
    chdir(glowpath)
    yd,utsec = datetime2yd(dt)[:2]

    #z = range(80,110+1,1)
    z = list(range(30,110+1,1))
    z += (
         [111.5,113.,114.5,116.] +
         list(chain(range(118,150+2,2),range(153,168+3,3),range(172,180+4,4),
                    range(185,205+5,5),range(211,223+6,6),range(230,244+7,7),
                    range(252,300+8,8),range(309,345+9,9),range(355,395+10,10),
                    range(406,428+11,11))) +
         [440,453,467,482,498,515,533,551] +
         list(range(570,950+20,20))
         )

    z = asarray(z)
#%% (1) setup flux at top of ionosphere
    ener,dE = glowfort.egrid()

    phitop = glowfort.maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)

    phi = hstack((ener[:,None],dE[:,None],phitop[:,None]))
#%% (2) msis,iri,glow model
    ion,ecalc,photI,ImpI,isr = glowfort.aurora(z,yd,utsec,glat,glon%360,
                                             f107a,f107,f107p,ap,phi)
#%% handle the outputs including common blocks
    zeta=glowfort.cglow.zeta.T #columns 11:20 are identically zero

    ver = DataFrame(index=z,
                    data=zeta[:,:11],
                    columns=[3371, 4278, 5200, 5577, 6300,7320,10400,3466,
                             7774, 8446,3726])
    photIon = DataFrame(index=z,
                   data=hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)),
                    columns=['photoIoniz','eImpactIoniz','ne',
                    'nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO'])

    isrparam = DataFrame(index=z,
                         data=isr,
                         columns=['ne','Te','Ti'])

    phitop = DataFrame(index=phi[:,0], #eV
                       data=phi[:,2],  #diffnumflux
                       columns=['diffnumflux'])
    zceta = glowfort.cglow.zceta.T

    return ver,photIon,isrparam,phitop,zceta
#%% plot
def plotaurora(phitop,ver,zceta,photIon,isr,dtime,glat,glon,E0):
#%% incident flux at top of ionosphere
    ax = figure().gca()
    ax.plot(phitop.index,phitop['diffnumflux'])
    ax.set_title('Incident Flux for $E_0={}$'.format(E0),fontsize='x-large')
    ax.set_xlabel('Beam Energy [eV]',fontsize='large')
    ax.set_ylabel('Flux',fontsize='large')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-4)
    ax.tick_params(axis='both',which='major',labelsize='medium')
#%% results of impacts
    fg,axs = subplots(1,4,sharey=True, figsize=(15,8))
    fg.suptitle('{} ({},{})'.format(dtime,glat,glon),fontsize='x-large')
    tight_layout(pad=3.2, w_pad=0.3)

    ax = axs[0]
    ax.plot(ver.values,ver.index)
    ax.set_xlabel('VER for $E_0={}$'.format(E0),fontsize='large')
    ax.set_ylabel('altitude [km]',fontsize='large')
    ax.set_ylim(top=400,bottom=ver.index[0])
    ax.set_xscale('log')
    ax.set_xlim(left=1e-4)
    ax.legend(ver.columns)
    ax.set_title('Volume emission rate',fontsize='x-large')

    ax = axs[1]
    ax.plot(photIon[['photoIoniz','eImpactIoniz']],photIon.index)
    ax.set_xlabel('ionization',fontsize='large')
    ax.set_xscale('log')
    ax.set_xlim(left=1e-1)
    ax.set_ylim(top=400)
    ax.legend(photIon.columns[:2])
    ax.set_title('Photo and e$^-$ impact ionization for $E_0={}$'.format(E0),fontsize='x-large')

    ax = axs[2]
    ax.semilogx(photIon[['ne','nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO']], photIon.index)
    ax.set_xlabel('Density',fontsize='large')
    ax.set_xscale('log')
    ax.set_xlim(left=1e-3)
    ax.set_ylim(top=400)
    ax.legend(photIon.columns[2:])
    ax.set_title('Electron and Ion Densities for $E_0={}$'.format(E0),fontsize='x-large')

    ax = axs[3]
    ax.semilogx(isr[['Te','Ti']], isr.index)
    ax.set_xlabel('Temperature [K]',fontsize='large')
    ax.legend(isr.columns[1:])
    ax.set_ylim(top=400)
    ax.set_title('Particle Temperature for $E_0={}$'.format(E0),fontsize='x-large')

    for a in axs:
        a.grid(True)
        a.tick_params(axis='both',which='major',labelsize='medium')

#%% total energy deposition vs. altitude
    fg,axs = subplots(1,2,sharey=True, figsize=(15,8))
    fg.suptitle('{} ({},{})'.format(dtime,glat,glon),fontsize='x-large')
    tight_layout(pad=3.2, w_pad=0.3)

    ax = axs[0]
    tez = glowfort.cglow.tez
    ax.plot(tez,ver.index)
    ax.set_xscale('log')
    ax.set_xlim(left=1e-1)
    ax.set_ylim(top=ver.index[-1],bottom=ver.index[0])
    ax.set_xlabel('Energy Deposited',fontsize='large')
    ax.set_ylabel('Altitude [km]',fontsize='large')
    ax.set_title('Total Energy Depostiion for $E_0={}$'.format(E0),fontsize='x-large')
#%% e^- impact ionization rates from ETRANS
    ax = axs[1]
    sion = glowfort.cglow.sion
    sion = DataFrame(index=ver.index,data=sion.T,columns=['O','O2','N2'])
    ax.plot(sion,ver.index)
    ax.set_xscale('log')
    ax.set_xlim(left=1e-6)
    ax.set_xlabel('e$^-$ impact ioniz. rate',fontsize='large')
    ax.set_title('electron impact ioniz. rates for $E_0={}$'.format(E0),fontsize='x-large')
    #ax.legend(True)
#%% constituants of per-wavelength VER
#    zcsum = zceta.sum(axis=-1)

    ax = figure().gca()
    for zc in rollaxis(zceta,1):
        ax.plot(ver.index,zc)
    ax.set_xlabel('emission constituants for $E_0={}$'.format(E0),fontsize='large')
    #ax.legend(True)
