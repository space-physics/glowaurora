from pathlib import Path
from numpy import rollaxis
from numpy.ma import masked_invalid
from pandas import DataFrame,Panel
from matplotlib.pyplot import figure, subplots,tight_layout,draw
from matplotlib.ticker import MultipleLocator #LogFormatterMathtext,
from matplotlib.colors import LogNorm

dymaj=50
dymin=10
dpi = 100

def _nicez(ax,zlim):
    ax.autoscale(True,axis='both',tight=True)
    ax.set_ylim(zlim)
    ax.yaxis.set_major_locator(MultipleLocator(dymaj))
    ax.yaxis.set_minor_locator(MultipleLocator(dymin))
    ax.grid(True,which='major',linewidth=1.)
    ax.grid(True,which='minor',linewidth=0.5)
    ax.tick_params(axis='both',which='major',labelsize='medium')

def plotaurora(phitop,ver,zceta,photIon,isr,sion,t,glat,glon,prate,lrate,tez,
               E0=None,flux=None,sza=None,zlim=(None,None),makeplot=None,odir=''):
    if makeplot is None:
        return

    if E0 and sza:
        titlend = '$E_0={:.0f}$ eV  SZA={:.1f}$^\circ$'.format(E0,sza)
    else:
        titlend = ''

    if isinstance(ver,Panel):
        z=ver.major_axis.values #[km], for convenience
    elif isinstance(ver,DataFrame):
        z=ver.index
#%% neutral background (MSIS) and Te,Ti (IRI-90)
    if not 'eig' in makeplot:
        fg,axs = subplots(1,2,sharey=True,figsize=(15,8))
        fg.suptitle('{} ({},{}) '.format(t,glat,glon)+ titlend)

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
    plotprodloss(phitop.index,z,prate['final'],lrate['final'],t,glat,glon,zlim,'',titlend)
#%% volume emission rate
    fg,axs = subplots(1,3,sharey=False, figsize=(15,8))
    fg.suptitle('{} ({},{}) '.format(t,glat,glon) + titlend)
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
        ind=[3371., 4278., 5200., 5577., 6300., 7320., 10400., 3466., 7774., 8446., 3726.]
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
        hi=a.pcolormesh(EKpcolor,z,masked_invalid(r.values),norm=LogNorm()) #pcolormesh canNOT handle nan at all!
        cb=fg.colorbar(hi,ax=a)
        cb.set_label('[cm$^{-3}$ s$^{-1}$ eV$^{-1}$]',labelpad=0)
        a.set_xscale('log')
        a.set_xlabel('Beam Energy [eV]')
#        a.legend(prod.minor_axis,loc='best')  # old, when plotting for each energy / input spectrum
        _nicez(a,zlim)

#%% NOTE: candidate for loading from gridaurora.plots instead
def writeplots(fg,plotprefix,E0,method,odir):
    odir = Path(odir)
    draw() #Must have this here or plot doesn't update in animation multiplot mode!
    #TIF was not faster and was 100 times the file size!
    #PGF is slow and big file,
    #RAW crashes
    #JPG no faster than PNG
    if 'png' in method:
        cn = (odir / (plotprefix + 'beam{:.0f}.png'.format(E0))).resolve()
        print('write {}'.format(cn))
        fg.savefig(cn,bbox_inches='tight',format='png',dpi=dpi)  # this is slow and async
