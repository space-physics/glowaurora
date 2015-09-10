#!/usr/bin/env python3
"""
Trivial example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
"""
from matplotlib.pyplot import figure, subplots,tight_layout
from pandas import DataFrame
from numpy import hstack,arange,append,array,rollaxis
from os import chdir
try:
    import seaborn
except:
    pass
#
#sys.path.append(os.path.dirname(os.path.abspath(__file__))) #enables fortran .dat files in site-packages module directory
import glowaurora
from glowaurora import glowfort
chdir(glowaurora.__path__[0])

def glowaurora(nbins,eflux,e0,iyd,utsec,glat,glon,f107a,f107,f107p,ap):
#%% temporarily use glow grid instead of our own
    ener,dE = glowfort.egrid()
    phitop = glowfort.maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)
    phi = hstack((ener[:,None],dE[:,None],phitop[:,None]))
#%% glow model
    z = arange(80,110+1,1)
    z = append(z,array([111.5,113.,114.5,116.,118.,120.,122.,124.,126., 128.,130.,132.,134.,136.,138.,140.,142.,144.,146., 148.,150.,153.,156.,159.,162.,165.,168.,172.,176., 180.,185.,190.,195.,200.,205.,211.,217.,223.,230.,237.,244.,252.,260.,268.,276.,284.,292.,300.,309., 318.,327.,336.,345.,355.,365.,375.,385.,395.,406., 417.,428.,440.,453.,467.,482.,498.,515.,533.,551., 570.,590.,610.,630.,650.,670.,690.,710.,730.,750., 770.,790.,810.,830.,850.,870.,890.,910.,930.,950.]))

    ion,ecalc,photI,ImpI,isr = glowfort.aurora(z,iyd,utsec,glat,glon%360,
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

    phitop = DataFrame(index=phi[:,0],
                       data=phi[:,2],
                       columns=['flux'])
    zceta = glowfort.cglow.zceta.T

    return ver,photIon,isrparam,phitop,zceta
#%% plot
def plotaurora(phitop,ver,zceta,photIon,isr,dtime,glat,glon):
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

#%% total energy deposition vs. altittude
    fg,axs = subplots(1,4,sharey=True, figsize=(15,8))
    fg.suptitle('{} ({},{})'.format(dtime,glat,glon),fontsize='x-large')
    tight_layout(pad=3.2, w_pad=0.3)

    ax = axs[0]
    tez = glowfort.cglow.tez
    ax.plot(tez,ver.index)
    ax.set_xlabel('Energy Deposited',fontsize='large')
    ax.set_ylabel('Altitude [km]',fontsize='large')
    ax.set_title('Total Energy Depostiion',fontsize='x-large')
#%% e^- impact ionization rates from ETRANS
    ax = axs[1]
    sion = glowfort.cglow.sion
    sion = DataFrame(index=ver.index,data=sion.T,columns=['O','O2','N2'])
    ax.plot(sion,ver.index)
    ax.set_xlabel('e$^-$ impact ioniz. rate',fontsize='large')
    ax.set_title('electron impact ioniz. rates',fontsize='x-large')
    #ax.legend(True)
#%% constituants of per-wavelength VER
#    zcsum = zceta.sum(axis=-1)

    ax = figure().gca()
    for zc in rollaxis(zceta,1):
        ax.plot(ver.index,zc)
    ax.set_xlabel('emission constituants',fontsize='large')
    #ax.legend(True)
