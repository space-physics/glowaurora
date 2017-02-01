from pathlib import Path
from numpy.ma import masked_invalid
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

#%% neutral background (MSIS) and Te,Ti (IRI-90)
    if not 'eig' in makeplot:
        fg,axs = subplots(1,2,sharey=True,figsize=(15,8))
        fg.suptitle('{} ({},{}) '.format(t,glat,glon)+ titlend)

        ind = ['nO','nO2','nN2','nNO']
        ax = axs[0]
        ax.semilogx(photIon.loc[:,ind], photIon.z_km)
        ax.set_xlabel('Number Density')
        ax.set_xscale('log')
        ax.set_xlim(left=1e1)
        ax.set_ylabel('Altitude [km]')
        _nicez(ax,zlim)
        ax.legend(ind)
        ax.set_title('Neutral Number Density')

        ind=['Te','Ti']
        ax = axs[1]
        ax.semilogx(isr.loc[:,ind], isr.z_km)
        ax.set_xlabel('Temperature [K]')
        ax.legend(ind)
        _nicez(ax,zlim)
        ax.set_xlim(100,10000)
        ax.set_title('Background Temperature')

        writeplots(fg,'bg_',E0,makeplot,odir)
#%% production and loss rates for species
    if 'eig' in makeplot:
        plotprodloss(prate.loc['final',...],
                     lrate.loc['final',...],t,glat,glon,zlim,'',titlend)
#%% volume emission rate
    fg,axs = subplots(1,3,sharey=False, figsize=(15,8))
    fg.suptitle('{} ({},{}) '.format(t,glat,glon) + titlend)
    tight_layout(pad=3.2, w_pad=0.6)

#%% incident flux at top of ionosphere
    ax = axs[0]
    ax.plot(phitop.eV,phitop,marker='.')

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
    ax.plot(ver.loc[...,ind].values.squeeze(), ver.z_km)
    ax.set_xlabel('Volume Emission Rate')
    ax.set_ylabel('altitude [km]')
    _nicez(ax,zlim)
    ax.set_xscale('log')
    if not 'eig' in makeplot:
        ax.set_xlim(1e-5,1e3)
    ax.legend(ind,loc='best')
    ax.set_title('Volume Emission Rate: Visible')
#%% ver invisible
    ind = [3371,7320,10400,3466,7774, 8446,3726]
    ax = axs[2]
    ax.plot(ver.loc[...,ind].values.squeeze(), ver.z_km)
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
        ax.semilogx(photIon.loc[:,ind], photIon.z_km)
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
        plotenerdep(tez,t,glat,glon,zlim,titlend)
#%% e^- impact ionization rates from ETRANS
        fg,axs = subplots(1,2,sharey=True, figsize=(15,8))
        fg.suptitle('{} ({},{})  '.format(t,glat,glon) + titlend)
        tight_layout(pad=3.2, w_pad=0.3)

        ind=['photoIoniz','eImpactIoniz']
        ax = axs[0]
        ax.plot(photIon.loc[:,ind],photIon.z_km)
        ax.set_xlabel('ionization')
        ax.set_xscale('log')
        ax.set_xlim(left=1e-1)
        _nicez(ax,zlim)
        ax.legend(ind)
        ax.set_title('Photo and e$^-$ impact ionization')
        ax.set_ylabel('Altitude [km]')

        ind=['O','O2','N2']
        ax = axs[1]
        ax.plot(sion.T, sion.z_km)
        ax.set_xscale('log')
        ax.set_xlim(1e-5,1e4)
        _nicez(ax,zlim)
        ax.set_xlabel('e$^-$ impact ioniz. rate')

        ax.set_title('electron impact ioniz. rates')
        ax.legend(ind)

        writeplots(fg,'ioniz_',E0,makeplot,odir)
#%% constituants of per-wavelength VER
#    zcsum = zceta.sum(axis=-1)
    if not 'eig' in makeplot:
        fg,axs = subplots(3,4,sharey=True,figsize=(15,8))
        for ax,zc,i in zip(axs.ravel(),
                           zceta.transpose('wavelength_nm','type','z_km'),
                           zceta.wavelength_nm):
            ax.plot(zc.T,zc.z_km)
            ax.set_xscale('log')
            #ax.set_xlabel('emission constituants  ' + titlend)
            ax.set_ylabel('Altitude [km]')
            ax.set_title('{} angstrom'.format(i.values))
            #ax.legend(True)

        writeplots(fg,'constit_',E0,makeplot,odir)

def plotenerdep(tez,t,glat,glon,zlim,titlend=''):
    fg= figure()
    ax = fg.gca()
    if tez.ndim==1:
        ax.plot(tez,tez.z_km)
        ax.set_xscale('log')
        ax.set_xlim(1e-1,1e6)
        ax.set_xlabel('Energy Deposited')
    else:
        ax.set_label('Beam Energy [eV]')
        hi=ax.pcolormesh(tez.eV,tez.z_km,masked_invalid(tez.values),norm=LogNorm())
        cb=fg.colorbar(hi,ax=ax)
        cb.set_label('Energy Deposited')

    _nicez(ax,zlim)
    ax.set_title('Total Energy Depostiion')
    ax.set_ylabel('altitude [km]')

def plotprodloss(prod,loss,t,glat,glon,zlim,titlbeg='',titlend=''):
    fg,ax = subplots(1,2,sharey=True,figsize=(15,8))
    fg.suptitle(titlbeg + ' Volume Production/Loss Rates   {} ({},{}) '.format(t,glat,glon)+ titlend)

    ax[0].set_title('Volume Production Rates')
    ax[0].set_ylabel('altitude [km]')

    ax[1].set_title('Volume Loss Rates')

    for a,R in zip(ax,[prod,loss]):
        try:
            hi=a.pcolormesh(R.eV.values,
                            R.z_km.values,
                            masked_invalid(R.values),
                            norm=LogNorm()) #pcolormesh canNOT handle nan at all!
            cb=fg.colorbar(hi,ax=a)
            cb.set_label('[cm$^{-3}$ s$^{-1}$ eV$^{-1}$]',labelpad=0)
            a.set_xscale('log')
            a.set_xlabel('Beam Energy [eV]')
    #        a.legend(prod.minor_axis,loc='best')  # old, when plotting for each energy / input spectrum
            _nicez(a,zlim)
        except TypeError as e:
            print('prodloss plot error    {}'.format(e))

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
