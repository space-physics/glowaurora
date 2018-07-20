#!/usr/bin/env python
"""
Example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
"""
from pathlib import Path
import logging
from datetime import datetime
import xarray
import numpy as np
from os import chdir
#
from sciencedates import datetime2yd, find_nearest
try:
    from gridaurora import getApF107
    from gridaurora.zglow import glowalt
except ImportError:
    getApF107 = glowalt = None
# Fortran GLOW
import glowfort
# path monkey patching
glowpath=Path(__file__).parents[1]
oldcwd = Path.cwd()

GRIDERR = 'gridaurora is required for this program.  pip install gridaurora.'


def runglowaurora(params:dict, z_km:np.ndarray=None) -> xarray.Dataset:
    """ Runs Fortran GLOW program and collects results in friendly arrays with metadata. """
    #%% (-2) check/process user inputs
    assert isinstance(params['flux'],(float,int,np.ndarray))
    assert isinstance(params['E0'],   (float,'float32',int))
    assert isinstance(params['t0'],   (datetime,str))
    assert isinstance(params['glat'], (float,int))
    assert isinstance(params['glon'], (float,int))
#%% (-1) if no manual f10.7 and ap, autoload by date
    if not 'f107a' in params or params['f107a'] is None:
        ifgetApF107 is None:
            raise ImportError(GRIDERR)
        f107Ap = getApF107(params['t0'])
        params['f107a'] = params['f107p'] = f107Ap['f107s'].item()
        params['f107']  = f107Ap['f107'].item()
        params['Ap']    = (f107Ap['Ap'].item(),)*7

    chdir(glowpath) #FIXME: hack for path issue
#%% flux grid / date
    eflux = np.atleast_1d(params['flux'])

    yeardoy,utsec = datetime2yd(params['t0'])[:2]
#%% (0) define altitude grid [km]
    if z_km is None:
        if glowalt is not None:
            z_km = glowalt()
        else:
            raise ImportError(GRIDERR)
#%% (1) setup flux at top of ionosphere
    ener,dE = glowfort.egrid()

    if eflux.size == 1:
        logging.info('generating maxwellian input differential number flux spectrum')
        # maxwellian input PhiTop at top of ionosphere
        phitop = glowfort.maxt(eflux,params['E0'],ener, dE, itail=0, fmono=0, emono=0)
    elif eflux.size > 1: #eigenprofile generation, one non-zero bin at a time
        logging.info('running in eigenprofile mode')
        e0ind = find_nearest(ener,params['e0'])[0] # FIXME should we interpolate instead? Maybe not, as long as we're consistent ref. Semeter 2006
        phitop = np.zeros_like(ener)
        phitop[e0ind] = ener[e0ind] #value in glow grid closest to zett grid
    else:
        return ValueError('I do not understand your electron flux input. Should be scalar or vector')

    phi = np.stack((ener,dE,phitop), 1)   # Nalt x 3
    assert phi.shape[1] == 3
#%% (2) msis,iri,glow model
    ion,ecalc,photI,ImpI,isr,prate,lrate,UV = glowfort.aurora(z_km,yeardoy,utsec,params['glat'],params['glon']%360,
                                                              params['f107a'],params['f107'],params['f107p'],params['Ap'],phi)

#%% (3) collect outputs
    zeta=glowfort.cglow.zeta.T #columns 11:20 are identically zero

    lamb=[3371., 4278., 5200., 5577., 6300., 7320., 10400., 3466., 7774., 8446., 3726.,
          1356., 1304., 1027., 989., 1900.]
    products=['nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+','nO','nO2','nN2','nNO']

    sim = xarray.Dataset()

    sim['ver'] = xarray.DataArray(
                   np.concatenate((zeta[:,:11],UV.T), axis=1),
                   dims=['z_km','wavelength_nm'],
                    coords={'z_km':z_km, 'wavelength_nm':lamb})

    sim['photIon'] = xarray.DataArray(np.hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)),
                               dims=['z_km','type'],
                               coords={'z_km':z_km,
                                'type':['photoIoniz','eImpactIoniz','ne']+products})

    sim['isr'] = xarray.DataArray(isr,
                          dims=['z_km','param'],
                          coords={'z_km':z_km,'param':['ne','Te','Ti']})

    sim['phitop'] = xarray.DataArray(phi[:,2],
                           dims=['eV'],
                           coords={'eV':phi[:,0]})

    sim.attrs['sza'] = np.degrees(glowfort.cglow.sza)

    sim['tez'] = xarray.DataArray(glowfort.cglow.tez,
                           dims=['z_km'], coords={'z_km':z_km})


#%% production and loss rates
    prate = prate.T; lrate=lrate.T #fortran to C ordering 2x170x20, only first 12 columns are used

    #column labels by inspection of fortran/gchem.f staring after "DO 150 I=1,JMAX" (thanks Stan!)
    sim['prates'] = xarray.DataArray(prate[1,:,:12], #columns 12:20 are identically zero
                      dims=['z_km','reaction'],
                      coords={'z_km':z_km,
                              'reaction':['O+(2P)','O+(2D)','O+(4S)','N+','N2+','O2+','NO+',
                                 'N2(A)','N(2P)','N(2D)','O(1S)','O(1D)']}
                    )

    sim['lrates'] = xarray.DataArray(lrate[1,:,:12], #columns 12:20 are identically zero
                      dims=['z_km','reaction'],
                      coords={'z_km':z_km,
                        'reaction':['O+(2P)','O+(2D)','O+(4S)','N+','N2+','O2+','NO+',
                                 'N2(A)','N(2P)','N(2D)','O(1S)','O(1D)']}
                    )

    sim['sion'] = xarray.DataArray(glowfort.cglow.sion,
                         dims=['gas','z_km'],
                         coords={'gas':['O','O2','N2'],
                                 'z_km':z_km})
# %% cleanup
    chdir(oldcwd)

    return sim


def rundayglow(t0,glat,glon,f107a,f107,f107p,ap, conj=True):
    '''
    Run GLOW for no auroral input, to simulate Dayglow.

    conj = whether to account for photoelectrons from conjugate hemisphere

    After running, extra variables can be found in glowfort.cglow. Unfortunately
    you have to dig into the Fortran code to see what they mean.
    '''

#%% (-2) check/process user inputs
    assert isinstance(t0,   (datetime,str))
    assert isinstance(glat, (float,int))
    assert isinstance(glon, (float,int))

    chdir(glowpath) #FIXME: hack for path issue

#%% flux grid / date

    yd,utsec = datetime2yd(t0)[:2]
#%% (0) define altitude grid [km]
    #z = glowalt()
    z = np.concatenate((range(30,110,1),np.logspace(np.log10(110),np.log10(1200),90)))

#%% (1) setup external flux at top of ionosphere. Set it to zero. Photoelectron flux from
#       conjugate hemisphere will be calculated internally, if conj=True.
    ener,dE = glowfort.egrid()
    phitop = np.zeros_like(ener)
    phi = np.hstack((ener[:,None],dE[:,None],phitop[:,None]))

#%% (2) msis,iri,glow model
    iconj = int(conj) # convert boolean to int for passing to Fortran.
    ion,ecalc,photI,ImpI,isr,UV = glowfort.dayglow(z,yd,utsec,glat,glon%360,
                                             f107a,f107,f107p,ap,phi,iconj)

#%% handle the outputs including common blocks
    zeta=glowfort.cglow.zeta.T #columns 11:20 are identically zero

    lamb=[3371., 4278., 5200., 5577., 6300., 7320., 10400., 3466., 7774., 8446., 3726.,
          1356., 1304., 1027., 989., 1900.]
    products=['nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+','nO','nO2','nN2','nNO']

    sim = xarray.Dataset()

    sim['ver'] = xarray.DataArray(np.concatenate((zeta[:,:11],UV.T), axis=1),
                   dims=['z_km','wavelength_nm'],
                    coords={'z_km':z,
                            'wavelength_nm':lamb})

    sim['photIon'] = xarray.DataArray(dims=['z_km','type'],
                        coords={'z_km':z,
                                'type':['photoIoniz','eImpactIoniz','ne']+products},
                   data=np.hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)))

    sim['isr'] = xarray.DataArray(dims=['z_km','param'],
                        coords={'z_km':z,'param':['ne','Te','Ti']},
                         data=isr)

    sim['phitop'] = xarray.DataArray(dims=['eV'],
                       coords={'eV':phi[:,0]},
                       data=phi[:,2])

    sim['zceta'] = xarray.DataArray(dims=['z_km','wavelength_nm','type'],
                      coords={'z_km':z, 'wavelength_nm':lamb[:11]},
                    data=glowfort.cglow.zceta.T[:,:11,:])  #Nalt x Nwavelengths  xNproductionEmissions


    sim.attrs['sza'] = np.degrees(glowfort.cglow.sza)

    sim['tez'] = xarray.DataArray(dims=['z_km'],
                    coords={'z_km':z},
                       data=glowfort.cglow.tez)

    sim['sion'] = xarray.DataArray(dims=['gas','z_km'],
                    coords={'gas':['O','O2','N2'],
                             'z_km':z},data=glowfort.cglow.sion)


    chdir(oldcwd)

    return sim


def verprodloss(params:dict):
    """ for a single time, computes VER, production, and loss vs. unit input flux
    inputs:
    -------
    t: a single datetime() when the eigenprofiles should be computed (solar zenith angle computed in Fortran code)
    glatlon: geographic coordinates of site (magnetic coordinates computed in Fortran code)
    flux: a vector of scaled differential number flux for each energy bin, scaled to "unit" energy flux for your eigenprofile input
    EK: a vector of energies [eV] to compute unit responses

    """
    simparams = ('prates','lrates','ver','tez')
    if isinstance(params['t0'],(list,tuple)):
        params['t0'] = params['t0'][0]
        print('using time',params['t0'])

    assert isinstance(params['t0'], datetime)

    EK = np.atleast_1d(params['EK'])

    sim = xarray.Dataset(coords={'EK':EK})

    for e in EK:
        print(f'{params["t0"]} E0: {e:.0f}')
        params['E0'] = e

        s = runglowaurora(params)

        if len(sim.variables)==1:
            for p in simparams:
                sim[p] = s[p].expand_dims('eV')
            sim.attrs['sza'] = s.sza
        else:
            for p in simparams:
                sim[p].loc[e,...] = s[p].expand_dims('eV')

        #plotaurora(sim)

    return sim


def ekpcolor(eigen):
    """ Loads electron energy bins from .csv file.
    This would be the gridding that a Transcar simulation used, so that Glow and Transcar
    are run on the same energy bins.
    """

    if isinstance(eigen, xarray.DataArray):
        e0 = eigen.loc[:,'low'].values
        eEnd = eigen.loc[:,'high'][-1]
        diffnumflux = eigen.loc[:,'flux'].values
    elif isinstance(eigen, (str,Path)):
        eigen = Path(eigen).expanduser()
        if eigen.suffix == '.csv':
            e0 =   np.loadtxt(eigen, usecols=[0], delimiter=',')
            eEnd = np.loadtxt(eigen, usecols=[1], delimiter=',')[-1]
            diffnumflux = None
        elif eigen.suffix == '.h5':
            import pandas
            bins = pandas.read_hdf(eigen)
            e0 = bins['low'].values
            eEnd = bins['high'].iloc[-1]
            diffnumflux = bins['flux'].values
        else:
            raise ValueError(f'I do not understand what file you want me to read {eigen}')
    else:
        raise ValueError(f'unknown data type {type(eigen)}')

    return np.append(e0,eEnd),e0,diffnumflux


def makeeigen(params:dict):
    """ Make eigenprofiles on an energy grid. Like Dahlgren 2012, Hirsch 2015, etc. """

    #for t in T:
    sim = verprodloss(params)

#TODO time stack

    return sim
