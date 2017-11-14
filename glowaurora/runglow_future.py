#!/usr/bin/env python
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
except ImportError:
    pass
#
from sciencedates import datetime2yd
from pyiri90.runiri90 import runiri
from msise00.runmsis import rungtd1d
import glowaurora
from glowaurora import glowfort
#
glowpath=glowaurora.__path__[0]

def runglowaurora(eflux,e0,dt,glat,glon,f107a,f107,f107p,ap,mass):
    chdir(glowpath)
    yd,utsec = datetime2yd(dt)[:2]

    z = arange(80,110+1,1)
    z = append(z,array([111.5,113.,114.5,116.,118.,120.,122.,124.,126., 128.,130.,132.,134.,136.,138.,140.,142.,144.,146., 148.,150.,153.,156.,159.,162.,165.,168.,172.,176., 180.,185.,190.,195.,200.,205.,211.,217.,223.,230.,237.,244.,252.,260.,268.,276.,284.,292.,300.,309., 318.,327.,336.,345.,355.,365.,375.,385.,395.,406., 417.,428.,440.,453.,467.,482.,498.,515.,533.,551., 570.,590.,610.,630.,650.,670.,690.,710.,730.,750., 770.,790.,810.,830.,850.,870.,890.,910.,930.,950.]))

    glowfort.cglow.zz = z*1e5
    glowfort.cglow.znd[:]=0.

# Set other parameters and switches:
    glowfort.cglow.jlocal = 0
    glowfort.cglow.kchem = 4
    glowfort.cglow.iscale = 1
    glowfort.cglow.xuvfac = 3.
    glowfort.cglow.hlybr = 0.
    glowfort.cglow.fexvir = 0.
    glowfort.cglow.hlya = 0.
    glowfort.cglow.heiew = 0.
#%% (1) setup flux at top of ionosphere
    ener,dE = glowfort.egrid()

    phitop = glowfort.maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)

    phi = hstack((ener[:,None],dE[:,None],phitop[:,None]))

    glowfort.cglow.zz = z*1e5
    glowfort.cglow.znd[:]=0.
#%% (2) call MSIS
    tselecopts = (1,)*25
    dens,temp=rungtd1d(dt,z,glat,glon,f107a,f107,ap,mass,tselecopts)
    glowfort.cglow.zo = dens['O']
    glowfort.cglow.zn2 = dens['N2']
    glowfort.cglow.zo2 = dens['O2']
    glowfort.cglow.zrho = dens['Total']
    glowfort.cglow.zns  = dens['N']
    glowfort.cglow.ztn  = temp['heretemp']
#%% (3) call snoem
    """
    Call SNOEMINT to obtain NO profile from the Nitric Oxide Empirical
    Model (NOEM)
    """
    glowfort.cglow.zno = glowfort.snoemint(dt.strftime('%Y%j'),
                               glat,glon,f107,ap,z,temp['heretemp'])
#%% (4a) call iri-90
    outf,oarr = runiri(dt,z,glat,glon,f107,f107a,ap,mass=48)
    chdir(glowpath) #need this since iri90 changes path
#%% (4b) store iri90 in COMMON blocks, after unit conversion
    glowfort.cglow.ze = outf['ne']/1e6 # M-3 -> CM-3
    glowfort.cglow.ze[glowfort.cglow.ze<100.] = 100.

    glowfort.cglow.zti = outf['Ti']
    i = glowfort.cglow.zti<glowfort.cglow.ztn
    glowfort.cglow.zti[i] = glowfort.cglow.ztn[i]

    glowfort.cglow.zte = outf['Te']
    i = glowfort.cglow.zte<glowfort.cglow.ztn
    glowfort.cglow.zte[i] = glowfort.cglow.ztn[i]

    glowfort.cglow.zxden[2,:] = outf['nO+']/1e6
    glowfort.cglow.zxden[5,:] = outf['nO2+']/1e6
    glowfort.cglow.zxden[6,:] = outf['nNO+']/1e6
#%% glow model

    ion,ecalc,photI,ImpI,isr = glowfort.aurora(z,yd,utsec,glat,glon%360,
                                             f107a,f107,phi)
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
def plotaurora(phitop,ver,zceta,photIon,isr,dtime,glat,glon):
#%% incident flux at top of ionosphere
    ax = figure().gca()
    ax.plot(phitop.index,phitop['diffnumflux'])
    ax.set_title('Incident Flux',fontsize='x-large')
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
    ax.set_xlabel('VER',fontsize='large')
    ax.set_ylabel('altitude [km]',fontsize='large')
    ax.set_ylim(top=ver.index[-1],bottom=ver.index[0])
    ax.set_xscale('log')
    ax.set_xlim(left=1e-4)
    ax.legend(ver.columns)
    ax.set_title('Volume emission rate',fontsize='x-large')

    ax = axs[1]
    ax.plot(photIon[['photoIoniz','eImpactIoniz']],photIon.index)
    ax.set_xlabel('ionization',fontsize='large')
    ax.set_xscale('log')
    ax.set_xlim(left=1e-1)
    ax.legend(photIon.columns[:2])
    ax.set_title('Photo and e$^-$ impact ionization',fontsize='x-large')

    ax = axs[2]
    ax.semilogx(photIon[['ne','nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO']], photIon.index)
    ax.set_xlabel('Density',fontsize='large')
    ax.set_xscale('log')
    ax.set_xlim(left=1e-3)
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
    ax.set_title('Total Energy Depostiion',fontsize='x-large')
#%% e^- impact ionization rates from ETRANS
    ax = axs[1]
    sion = glowfort.cglow.sion
    sion = DataFrame(index=ver.index,data=sion.T,columns=['O','O2','N2'])
    ax.plot(sion,ver.index)
    ax.set_xscale('log')
    ax.set_xlim(left=1e-6)
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
