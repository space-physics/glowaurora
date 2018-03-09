from pathlib import Path
import xarray
import logging
from numpy.ma import masked_invalid
from matplotlib.pyplot import figure, draw
from matplotlib.ticker import MultipleLocator #LogFormatterMathtext,
from matplotlib.colors import LogNorm

dymaj=50
dymin=10
dpi = 100

def _nicez(ax, params:dict):
    """ make ordinate axis look nice for altitude """
    ax.autoscale(True,axis='both',tight=True)
    ax.set_ylim(params['zlim'])
    ax.yaxis.set_major_locator(MultipleLocator(dymaj))
    ax.yaxis.set_minor_locator(MultipleLocator(dymin))
    ax.grid(True,which='major',linewidth=1.)
    ax.grid(True,which='minor',linewidth=0.5)
    ax.tick_params(axis='both',which='major',labelsize='medium')

def _plotisr(params,sim,st:str):

    if not 'photIon' in sim or not 'isr' in sim:
        return

    fg = figure(figsize=(15,8))
    axs = fg.subplots(1,2,sharey=True,)
    fg.suptitle(st)

    ind = ['nO','nO2','nN2','nNO']
    ax = axs[0]

    ax.semilogx(sim['photIon'].loc[:,ind],
                sim.z_km)
    ax.set_xlabel('Number Density')
    ax.set_xscale('log')
    ax.set_xlim(left=1e1)
    ax.set_ylabel('Altitude [km]')
    _nicez(ax,params)
    ax.legend(ind)
    ax.set_title('Neutral Number Density')

    ind=['Te','Ti']
    ax = axs[1]
    # Need to use .values here.
    ax.semilogx(sim['isr'].loc[:,ind],
                sim.z_km)
    ax.set_xlabel('Temperature [K]')
    ax.legend(ind)
    _nicez(ax,params)
    ax.set_xlim(100,10000)
    ax.set_title('Background Temperature')

    writeplots(fg,'bg_',params)

def _plotver(sim,params,supertitle):
    if not 'phitop' in sim:
        return

    fg = figure(figsize=(15,8))
    axs = fg.subplots(1,3,sharey=False)
    fg.suptitle(supertitle)
    fg.tight_layout(pad=3.2, w_pad=0.6)

#%% incident flux at top of ionosphere
    ax = axs[0]
    ax.plot(sim['phitop'].eV, sim['phitop'], marker='.')

    titxt=f'Total Incident Flux={params["flux"]:.1f},'
    ax.set_title(titxt)

    ax.set_xlabel('Beam Energy [eV]')
    ax.set_ylabel('Flux [erg sr$^{-1}$ s$^{-1}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-4,1e6)
    ax.grid(True)
# %% ver visible
    if not 'ver' in sim:
        return

    ind= [4278, 5200, 5577, 6300]
    ax = axs[1]
    ax.plot(sim['ver'].loc[...,ind],
            sim.z_km)
    ax.set_xlabel('Volume Emission Rate')
    ax.set_ylabel('altitude [km]')
    _nicez(ax,params)
    ax.set_xscale('log')
    if not 'eig' in params['makeplot']:
        ax.set_xlim(1e-5,1e3)
    ax.legend(ind,loc='best')
    ax.set_title('Volume Emission Rate: Visible')
#%% ver invisible
    ind = [3371,7320,10400,3466,7774, 8446,3726,1356., 1304., 1027., 989., 1900.]
    ax = axs[2]
    ax.plot(sim['ver'].loc[...,ind],
            sim.z_km)
    ax.set_xlabel('Volume Emission Rate')
    _nicez(ax,params)
    ax.set_xscale('log')
    if not 'eig' in params['makeplot']:
        ax.set_xlim(1e-5,1e3)
    ax.legend(ind,loc='best')
    ax.set_title('Volume Emission Rate: IR & UV')

    writeplots(fg,'ver_',params)


def _plotdens(sim,params,supertitle):

    if not 'photIon' in sim:
        return

    ind=['ne','nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+']
    fg=figure()
    ax = fg.gca()
    ax.semilogx(sim['photIon'].loc[:,ind],
                sim.z_km)
    ax.set_xlabel('Density')
    ax.set_xscale('log')
    ax.set_xlim(left=1e-3)
    _nicez(ax,params)
    ax.legend(ind)
    ax.set_title('Electron and Ion Densities')
    ax.set_ylabel('Altitude [km]')

    writeplots(fg,'effects_',params)


def _plotioniz(sim,params,supertitle):
    if 'photIon' not in sim:
        return

    fg = figure(figsize=(15,8))
    axs = fg.subplots(1,2,sharey=True,)
    fg.suptitle(supertitle)
    fg.tight_layout(pad=3.2, w_pad=0.3)

    ind=['photoIoniz','eImpactIoniz']
    ax = axs[0]
    ax.plot(sim['photIon'].loc[:,ind],
            sim.z_km)
    ax.set_xlabel('ionization')
    ax.set_xscale('log')
    ax.set_xlim(left=1e-1)
    _nicez(ax,params)
    ax.legend(ind)
    ax.set_title('Photo and e$^-$ impact ionization')
    ax.set_ylabel('Altitude [km]')

    ind=['O','O2','N2']
    ax = axs[1]
    ax.plot(sim['sion'].T,
            sim.z_km)
    ax.set_xscale('log')
    ax.set_xlim(1e-5,1e4)
    _nicez(ax,params)
    ax.set_xlabel('e$^-$ impact ioniz. rate')

    ax.set_title('electron impact ioniz. rates')
    ax.legend(ind)

    writeplots(fg,'ioniz_',params)


def plotaurora(params:dict, sim:xarray.Dataset):
    """ Plot all sorts of auroral/dayglow parameters from GLOW simulation. """
    makeplot = params['makeplot']
    if not makeplot:
        return

    if 'EK' in params:
        EK = params['EK']
    else:
        EK = params['E0']

    titlend = f'$E_0={EK}$ eV  SZA={sim.sza:.1f}$^\circ$'
    supertitle = f'{params["t0"]} ({params["glat"]},{params["glon"]}) \n {titlend}'

#%% neutral background (MSIS) and Te,Ti (IRI-90)
    if not 'eig' in makeplot:
        _plotisr(params,sim,supertitle)
#%% production and loss rates for species
    plotprodloss(sim['prates'], sim['lrates'], params, supertitle)
#%% volume emission rate
    _plotver(sim,params,supertitle)
# %% Ne, Ni
    if not 'eig' in makeplot:
        _plotdens(sim,params,supertitle)
# %% total energy deposition vs. altitude
    if not 'eig' in makeplot:
        plotenerdep(sim['tez'],params,titlend)
#% % e^- impact ionization rates from ETRANS
        _plotioniz(sim,params,supertitle)
#%% constituants of per-wavelength VER
    if not 'eig' in makeplot:
        _plotconstit(sim,params,supertitle)
\
def _plotconstit(sim,params,supertitle):
    # FIXME: talk to Stan about what these  mean, etc.
    if 'zceta' not in sim:
        return

    zcsum = sim['zceta'].sum(axis=-1)

    fg = figure(figsize=(15,8))
    axs = fg.subplots(3,4,sharey=True)
    for ax,zc,i in zip(axs.ravel(),
                       sim['zceta'].transpose('wavelength_nm','type','z_km'),
                       sim.wavelength_nm):
        ax.plot(zc.T,zc.z_km)
        ax.set_xscale('log')
        #ax.set_xlabel('emission constituants  ' + titlend)
        ax.set_ylabel('Altitude [km]')
        ax.set_title(f'{i.values} angstrom')
        #ax.legend(True)

    writeplots(fg,'constit_',params)


def plotenerdep(tez,params,titlend=''):
    """ plot energy deposition vs. altitude """
    fg= figure()
    ax = fg.gca()

    tez = tez.squeeze()

    if tez.ndim==1:
        ax.plot(tez, tez.z_km)
        ax.set_xscale('log')
        ax.set_xlim(1e-1,1e6)
        ax.set_xlabel('Energy Deposited')
    else:
        ax.set_label('Beam Energy [eV]')
        hi=ax.pcolormesh(tez.eV,tez.z_km,masked_invalid(tez.values),norm=LogNorm())
        cb=fg.colorbar(hi,ax=ax)
        cb.set_label('Energy Deposited')

    _nicez(ax,params)
    ax.set_title('Total Energy Depostiion')
    ax.set_ylabel('altitude [km]')


def plotprodloss(prod,loss,params,st):
    """ plot production/loss vs. alttiude """

    fg = figure(figsize=(15,8))
    ax = fg.subplots(1,2,sharey=True)
    fg.suptitle(f' Volume Production/Loss Rates   {st}')

    for a,R,title in zip(ax,[prod,loss],('Production','Loss')):
        if R.ndim==2:
            a.set_title('Volume {title} Rates')
            try:
                hi=a.pcolormesh(R.values,
                                R.z_km.values,
                                masked_invalid(R.values),
                                norm=LogNorm()) #pcolormesh canNOT handle nan at all!
                cb=fg.colorbar(hi,ax=a)
                cb.set_label('[cm$^{-3}$ s$^{-1}$ eV$^{-1}$]',labelpad=0)
                a.set_xscale('log')
                a.set_xlabel('Beam Energy [eV]')
                ax[0].set_ylabel('altitude [km]')
                _nicez(a,params)
            except TypeError as e:
                logging.warning(f'prodloss plot error    {e}')
        elif R.ndim==3:
            a.set_title('Volume {title} Rates')
            for a,R in zip(ax,[prod,loss]):
                for r in R:
                    for s in r.T:
                        a.plot(s, s.z_km, label=s.reaction.item())

                a.set_xlabel('[cm$^{-3}$ s$^{-1}$ eV$^{-1}$]',labelpad=0)
                a.set_ylabel('altitude [km]')
                a.set_xscale('log')
                a.legend(loc='best')
                _nicez(a,params)

# %%
def writeplots(fg:figure, plotprefix:str, params:dict):
    """ Save Matplotlib plots to disk.

    inputs:
    ------

    fg: Matplotlib figure handle
    odir: directory to write image into e.g. for this particular simulation.
    plotprefix: stem of filename
    E0: beam energy (eV)
    method: format of image

    Some observations on image formats:

      * TIF was not faster and was 100 times the file size!
      * PGF is slow and big file,
      * RAW crashes
      * JPG no faster than PNG
    """
    if 'odir' not in params or not params['odir']:
        return

    odir = Path(params['odir']).expanduser()
    draw() #Must have this here or plot doesn't update in animation multiplot mode!

    cn = odir / (plotprefix + f'beam{params["E0"]:.0f}.{params["plotformat"]}')
    print('write',cn)
    fg.savefig(cn, bbox_inches='tight', dpi=dpi)  # this is slow and async
