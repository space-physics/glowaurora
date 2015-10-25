#!/usr/bin/env python3
"""
Example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
"""
from __future__ import division,absolute_import
import logging
from datetime import datetime
from six import integer_types,string_types
from itertools import chain
from matplotlib.pyplot import figure, subplots,tight_layout,draw
from matplotlib.ticker import MultipleLocator #LogFormatterMathtext,
from matplotlib.colors import LogNorm
from pandas import DataFrame,Panel
from numpy import hstack,asarray,rollaxis,degrees,zeros_like,ndarray,atleast_1d
from numpy.ma import masked_invalid
from os import chdir
from os.path import join
try:
    import seaborn
    seaborn.set_set_context('poster')
except:
    pass
#
from histutils.fortrandates import datetime2yd
from histutils.findnearest import find_nearest
import glowaurora
from glowaurora import glowfort
#
glowpath=glowaurora.__path__[0]

dymaj=50
dymin=10
dpi = 100

def runglowaurora(eflux,e0,t0,glat,glon,f107a,f107,f107p,ap):
#%% (-1) check/process user inputs
    assert isinstance(eflux,(float,integer_types,ndarray))
    assert isinstance(e0,   (float,integer_types))
    assert isinstance(t0,   (datetime,string_types))
    assert isinstance(glat, (float,integer_types))
    assert isinstance(glon, (float,integer_types))
    assert isinstance(f107a,(float,integer_types))
    assert isinstance(f107, (float,integer_types))
    assert isinstance(f107p,(float,integer_types))
    assert isinstance(ap,   (float,integer_types))

    eflux = atleast_1d(eflux)
    chdir(glowpath) #FIXME: hack for path issue
    yd,utsec = datetime2yd(t0)[:2]
#%% (0) define altitude grid [km] (here i liberally copied that used by Stan)
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

    if len(eflux)==1:
        logging.info('generating maxwellian input differential number flux spectrum')
        # maxwellian input PhiTop at top of ionosphere
        phitop = glowfort.maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)
    elif len(eflux)>1: #eigenprofile generation, one non-zero bin at a time
        logging.info('running in eigenprofile mode')
        e0ind = find_nearest(ener,e0)[0] #FIXME should we interpolate instead? Maybe not, as long as we're consistent ref. Semeter 2006
        phitop = zeros_like(ener)
        phitop[e0ind] = ener[e0ind] #value in glow grid closest to zett grid
    else:
        return TypeError('I do not understand your electron flux input. Should be scalar or vector')

    phi = hstack((ener[:,None],dE[:,None],phitop[:,None]))
#%% (2) msis,iri,glow model
    ion,ecalc,photI,ImpI,isr,prate,lrate = glowfort.aurora(z,yd,utsec,glat,glon%360,
                                             f107a,f107,f107p,ap,phi)
#%% handle the outputs including common blocks
    zeta=glowfort.cglow.zeta.T #columns 11:20 are identically zero

    ver = DataFrame(index=z,
                    data=zeta[:,:11],
                    columns=[3371, 4278, 5200, 5577, 6300,7320,10400,3466,7774, 8446,3726])
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

    zceta = glowfort.cglow.zceta.T  #Nalt x Nwavelengths  xNproductionEmissions

    sza = degrees(glowfort.cglow.sza)

    tez = glowfort.cglow.tez #1-D vs. altitude

#%% production and loss rates
    prate = prate.T; lrate=lrate.T #fortran to C ordering 2x170x20, only first 12 columns are used

    #column labels by inspection of fortran/gchem.f staring after "DO 150 I=1,JMAX" (thanks Stan!)
    prates = Panel(items=['pre','final'],
                    major_axis=z,
                    minor_axis=['O+(2P)','O+(2D)','O+(4S)','N+','N2+','O2+','NO+',
                                 'N2(A)','N(2P)','N(2D)','O(1S)','O(1D)'],
                    data=prate[...,:12], #columns 12:20 are identically zero
                        )

    lrates = Panel(items=['pre','final'],
                    major_axis=z,
                    minor_axis=['O+(2P)','O+(2D)','O+(4S)','N+','N2+','O2+','NO+',
                                 'N2(A)','N(2P)','N(2D)','O(1S)','O(1D)'],
                    data=lrate[...,:12], #columns 12:20 are identically zero
                        )


    return ver,photIon,isrparam,phitop,zceta,sza,prates,lrates,tez

#%% plot
def _nicez(ax,zlim):
    ax.autoscale(True,axis='both',tight=True)
    ax.set_ylim(zlim)
    ax.yaxis.set_major_locator(MultipleLocator(dymaj))
    ax.yaxis.set_minor_locator(MultipleLocator(dymin))
    ax.grid(True,which='major',linewidth=1.)
    ax.grid(True,which='minor',linewidth=0.5)
    ax.tick_params(axis='both',which='major',labelsize='medium')

def plotaurora(phitop,ver,zceta,photIon,isr,dtime,glat,glon,prate,lrate,tez,
               E0=None,flux=None,sza=None,zlim=(None,None),makeplot=None,odir=''):
    if makeplot is None:
        return

    if E0 and sza:
        titlend = '$E_0={:.0f}$ eV  SZA={:.1f}$^\circ$'.format(E0,sza)
    else:
        titlend = ''

    z=ver.major_axis.values #[km], for convenience
#%% neutral background (MSIS) and Te,Ti (IRI-90)
    if not 'eig' in makeplot:
        fg,axs = subplots(1,2,sharey=True,figsize=(15,8))
        fg.suptitle('{} ({},{}) '.format(dtime,glat,glon)+ titlend)

        ind = ['nO','nO2','nN2','nNO']
        ax = axs[0]
        ax.semilogx(photIon[ind], z)
        ax.set_xlabel('Number Density')
        ax.set_xscale('log')
        ax.set_xlim(left=1e1)
        ax.set_ylabel('Altitude [km]')
        _nicez(ax,zlim)
        ax.legend(ind)
        ax.set_title('Neutral Number Density')

        ind=['Te','Ti']
        ax = axs[1]
        ax.semilogx(isr[ind], z)
        ax.set_xlabel('Temperature [K]')
        ax.legend(ind)
        _nicez(ax,zlim)
        ax.set_title('Background Temperature')

        writeplots(fg,'bg_',E0,makeplot,odir)
#%% production and loss rates for species
    plotprodloss(z,prate,lrate,dtime,glat,glon,zlim,'',titlend)
#%% volume emission rate
    fg,axs = subplots(1,3,sharey=False, figsize=(15,8))
    fg.suptitle('{} ({},{}) '.format(dtime,glat,glon) + titlend)
    tight_layout(pad=3.2, w_pad=0.6)

# incident flux at top of ionosphere
    ax = axs[0]
    ax.plot(phitop.index,phitop['diffnumflux'],marker='.')

    titxt='Incident Flux'
    if flux:
        titxt+='  Total Flux={:.1f},'.format(flux)
    ax.set_title(titxt)

    ax.set_xlabel('Beam Energy [eV]')
    ax.set_ylabel('Flux [erg sr$^{-1}$ s$^{-1}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-4,1e6)
    ax.grid(True)
# ver visible
    ind= [4278, 5200, 5577, 6300]
    ax = axs[1]
    ax.plot(ver.loc[:,:,ind].values.squeeze(),z)
    ax.set_xlabel('Volume Emission Rate')
    ax.set_ylabel('altitude [km]')
    _nicez(ax,zlim)
    ax.set_xscale('log')
    if not 'eig' in makeplot:
        ax.set_xlim(1e-5,1e3)
    ax.legend(ind,loc='best')
    ax.set_title('Volume Emission Rate: Visible')
# ver invisible
    ind = [3371,7320,10400,3466,7774, 8446,3726]
    ax = axs[2]
    ax.plot(ver.loc[:,:,ind].values.squeeze(),z)
    ax.set_xlabel('Volume Emission Rate')
    _nicez(ax,zlim)
    ax.set_xscale('log')
    if not 'eig' in makeplot:
        ax.set_xlim(1e-5,1e3)
    ax.legend(ind,loc='best')
    ax.set_title('Volume Emission Rate: IR & UV')

    writeplots(fg,'ver_',E0,makeplot,odir)
#%% Ne, Ni
    if not 'eig' in makeplot:
        ind=['ne','nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+']
        fg=figure()
        ax = fg.gca()
        ax.semilogx(photIon[ind], photIon.index)
        ax.set_xlabel('Density')
        ax.set_xscale('log')
        ax.set_xlim(left=1e-3)
        _nicez(ax,zlim)
        ax.legend(ind)
        ax.set_title('Electron and Ion Densities')
        ax.set_ylabel('Altitude [km]')

        writeplots(fg,'effects_',E0,makeplot,odir)
#%% total energy deposition vs. altitude
    if not 'eig' in makeplot:
        plotenerdep(z,tez,t,glat,glon,zlim,titlend)
#%% e^- impact ionization rates from ETRANS
        fg,axs = subplots(1,2,sharey=True, figsize=(15,8))
        fg.suptitle('{} ({},{})  '.format(t,glat,glon) + titlend)
        tight_layout(pad=3.2, w_pad=0.3)

        ind=['photoIoniz','eImpactIoniz']
        ax = axs[0]
        ax.plot(photIon[ind],photIon.index)
        ax.set_xlabel('ionization')
        ax.set_xscale('log')
        ax.set_xlim(left=1e-1)
        _nicez(ax,zlim)
        ax.legend(ind)
        ax.set_title('Photo and e$^-$ impact ionization')

        ind=['O','O2','N2']
        ax = axs[1]
        sion = glowfort.cglow.sion
        sion = DataFrame(index=z,data=sion.T,columns=ind)
        ax.plot(sion,z)
        ax.set_xscale('log')
        ax.set_xlim(1e-5,1e4)
        _nicez(ax,zlim)
        ax.set_xlabel('e$^-$ impact ioniz. rate')
        ax.set_ylabel('Altitude [km]')
        ax.set_title('electron impact ioniz. rates')
        ax.legend(ind)

        writeplots(fg,'ioniz_',E0,makeplot,odir)
#%% constituants of per-wavelength VER
#    zcsum = zceta.sum(axis=-1)
    if not 'eig' in makeplot:
        ind=[3371, 4278, 5200, 5577, 6300, 7320, 10400, 3466, 7774, 8446, 3726]
        fg,axs = subplots(3,4,sharey=True,figsize=(15,8))
        for ax,zc,i in zip(axs.ravel(),rollaxis(zceta,1)[:11,...],ind):
            ax.plot(zc,z)
            ax.set_xscale('log')
            #ax.set_xlabel('emission constituants  ' + titlend)
            ax.set_ylabel('Altitude [km]')
            ax.set_title('{} angstrom'.format(i))
            #ax.legend(True)

        writeplots(fg,'constit_',E0,makeplot,odir)

def plotenerdep(EKpcolor,z,tez,t,glat,glon,zlim,titlend=''):
    fg= figure()
    ax = fg.gca()
    if tez.ndim==1:
        ax.plot(tez,z)
        ax.set_xscale('log')
        ax.set_xlim(1e-1,1e6)
        ax.set_xlabel('Energy Deposited')
    else:
        ax.set_label('Beam Energy [eV]')
        hi=ax.pcolormesh(EKpcolor,z,masked_invalid(tez),norm=LogNorm())
        cb=fg.colorbar(hi,ax=ax)
        cb.set_label('Energy Deposited')

    _nicez(ax,zlim)
    ax.set_title('Total Energy Depostiion')

def plotprodloss(EKpcolor,z,prod,loss,t,glat,glon,zlim,titlbeg='',titlend=''):
    fg,ax = subplots(1,2,sharey=True,figsize=(15,8))
    fg.suptitle(titlbeg + ' Volume Production/Loss Rates   {} ({},{}) '.format(t,glat,glon)+ titlend)

    ax[0].set_title('Volume Production Rates')
    ax[0].set_ylabel('altitude [km]')

    ax[1].set_title('Volume Loss Rates')

    for a,r in zip(ax,[prod,loss]):
        hi=a.pcolormesh(EKpcolor,z,masked_invalid(r),norm=LogNorm()) #pcolormesh canNOT handle nan at all!
        cb=fg.colorbar(hi,ax=a)
        cb.set_label('[cm$^{-3}$ s$^{-1}$ eV$^{-1}$]',labelpad=0)
        a.set_xscale('log')
        a.set_xlabel('Beam Energy [eV]')
#        a.legend(prod.minor_axis,loc='best')  # old, when plotting for each energy / input spectrum
        _nicez(a,zlim)

#%%
def writeplots(fg,plotprefix,E0,method,odir):
    draw() #Must have this here or plot doesn't update in animation multiplot mode!
    #TIF was not faster and was 100 times the file size!
    #PGF is slow and big file,
    #RAW crashes
    #JPG no faster than PNG
    if 'png' in method:
        cn = join(odir,(plotprefix + 'beam{:.0f}.png'.format(E0)))
        print('write {}'.format(cn))
        fg.savefig(cn,bbox_inches='tight',format='png',dpi=dpi)  # this is slow and async
