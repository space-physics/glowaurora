#!/usr/bin/env python
# -*- coding: future_fstrings -*-
"""
Example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
"""
from pathlib import Path
import logging
from datetime import datetime
from xarray import DataArray
import numpy as np
from pandas import read_hdf
from os import chdir
#
from sciencedates import datetime2yd, find_nearest
from gridaurora import readmonthlyApF107
from gridaurora.zglow import glowalt
import glowfort
glowpath=Path(__file__).parents[0]
oldcwd = Path.cwd()

def runglowaurora(eflux,e0,t0,glat,glon,
                  f107a=None,f107=None,f107p=None,ap=None,f107apfn=None):


#%% (-2) check/process user inputs
    assert isinstance(eflux,(float,int,np.ndarray))
    assert isinstance(e0,   (float,'float32',int))
    assert isinstance(t0,   (datetime,str))
    assert isinstance(glat, (float,int))
    assert isinstance(glon, (float,int))
#%% (-1) if no manual f10.7 and ap, autoload by date
    if not(f107a and f107 and f107p and ap):
        f107Ap=readmonthlyApF107(t0,f107apfn)
        f107a = f107p = f107Ap['f107s']
        f107  = f107Ap['f107o']
        ap    = (f107Ap['Apo'],)*7

    chdir(str(glowpath)) #FIXME: hack for path issue

#%% flux grid / date

    eflux = np.atleast_1d(eflux)

    yd,utsec = datetime2yd(t0)[:2]
#%% (0) define altitude grid [km]
    z=glowalt()
#%% (1) setup flux at top of ionosphere
    ener,dE = glowfort.egrid()

    if len(eflux)==1:
        logging.info('generating maxwellian input differential number flux spectrum')
        # maxwellian input PhiTop at top of ionosphere
        phitop = glowfort.maxt(eflux,e0,ener, dE, itail=0, fmono=0, emono=0)
    elif len(eflux)>1: #eigenprofile generation, one non-zero bin at a time
        logging.info('running in eigenprofile mode')
        e0ind = find_nearest(ener,e0)[0] #FIXME should we interpolate instead? Maybe not, as long as we're consistent ref. Semeter 2006
        phitop = np.zeros_like(ener)
        phitop[e0ind] = ener[e0ind] #value in glow grid closest to zett grid
    else:
        return ValueError('I do not understand your electron flux input. Should be scalar or vector')

    phi = np.hstack((ener[:,None],dE[:,None],phitop[:,None]))
#%% (2) msis,iri,glow model
    ion,ecalc,photI,ImpI,isr,prate,lrate,UV = glowfort.aurora(z,yd,utsec,glat,glon%360,
                                             f107a,f107,f107p,ap,phi)

#%% handle the outputs including common blocks
    zeta=glowfort.cglow.zeta.T #columns 11:20 are identically zero

    lamb=[3371., 4278., 5200., 5577., 6300., 7320., 10400., 3466., 7774., 8446., 3726.,
          1356., 1304., 1027., 989., 1900.]
    products=['nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+','nO','nO2','nN2','nNO']

    ver = DataArray(dims=['z_km','wavelength_nm'],
                    coords={'z_km':z,
                            'wavelength_nm':lamb},
                    data=np.concatenate((zeta[:,:11],UV.T), axis=1))
    
    photIon = DataArray(dims=['z_km','type'],
                        coords={'z_km':z,
                                'type':['photoIoniz','eImpactIoniz','ne']+products},
                   data=np.hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)))

    isrparam = DataArray(dims=['z_km','param'],
                        coords={'z_km':z,'param':['ne','Te','Ti']},
                         data=isr)

    phitop = DataArray(dims=['eV'],
                       coords={'eV':phi[:,0]},
                       data=phi[:,2])

    zceta = DataArray(dims=['z_km','wavelength_nm','type'],
                      coords={'z_km':z, 'wavelength_nm':lamb[:11]},
                    data=glowfort.cglow.zceta.T[:,:11,:])  #Nalt x Nwavelengths  xNproductionEmissions

    # FIXME what is 3rd axis?


    sza = np.degrees(glowfort.cglow.sza)

    tez = DataArray(dims=['z_km'],
                    coords={'z_km':z},
                       data=glowfort.cglow.tez)


#%% production and loss rates
    prate = prate.T; lrate=lrate.T #fortran to C ordering 2x170x20, only first 12 columns are used

    #column labels by inspection of fortran/gchem.f staring after "DO 150 I=1,JMAX" (thanks Stan!)
    prates = DataArray(data=prate[...,:12], #columns 12:20 are identically zero
                      dims=['type','z_km','reaction'],
                      coords={'type':['pre','final'],
                               'z_km':z,
                        'reaction':['O+(2P)','O+(2D)','O+(4S)','N+','N2+','O2+','NO+',
                                 'N2(A)','N(2P)','N(2D)','O(1S)','O(1D)']}
                    )

    lrates = DataArray(data=lrate[...,:12], #columns 12:20 are identically zero
                      dims=['type','z_km','reaction'],
                      coords={'type':['pre','final'],
                               'z_km':z,
                        'reaction':['O+(2P)','O+(2D)','O+(4S)','N+','N2+','O2+','NO+',
                                 'N2(A)','N(2P)','N(2D)','O(1S)','O(1D)']}
                    )

    sion = DataArray(dims=['gas','z_km'],
                    coords={'gas':['O','O2','N2'],
                             'z_km':z},data=glowfort.cglow.sion)


    chdir(str(oldcwd))

    return ver,photIon,isrparam,phitop,zceta,sza,prates,lrates,tez,sion





def rundayglow(t0,glat,glon,f107a,f107,f107p,ap,conj=True):
    '''
    Run GLOW for no auroral input.
    conj = whether to account for photoelectrons from conjugate hemisphere
    
    After running, extra variables can be found in glowfort.cglow. Unfortunately
    you have to dig into the Fortran code to see what they mean.
    '''
    
#%% (-2) check/process user inputs
    assert isinstance(t0,   (datetime,str))
    assert isinstance(glat, (float,int))
    assert isinstance(glon, (float,int))

    chdir(str(glowpath)) #FIXME: hack for path issue

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

    ver = DataArray(dims=['z_km','wavelength_nm'],
                    coords={'z_km':z,
                            'wavelength_nm':lamb},
                    data=np.concatenate((zeta[:,:11],UV.T), axis=1))
    
    photIon = DataArray(dims=['z_km','type'],
                        coords={'z_km':z,
                                'type':['photoIoniz','eImpactIoniz','ne']+products},
                   data=np.hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)))

    isrparam = DataArray(dims=['z_km','param'],
                        coords={'z_km':z,'param':['ne','Te','Ti']},
                         data=isr)

    phitop = DataArray(dims=['eV'],
                       coords={'eV':phi[:,0]},
                       data=phi[:,2])

    zceta = DataArray(dims=['z_km','wavelength_nm','type'],
                      coords={'z_km':z, 'wavelength_nm':lamb[:11]},
                    data=glowfort.cglow.zceta.T[:,:11,:])  #Nalt x Nwavelengths  xNproductionEmissions


    sza = np.degrees(glowfort.cglow.sza)

    tez = DataArray(dims=['z_km'],
                    coords={'z_km':z},
                       data=glowfort.cglow.tez)

    sion = DataArray(dims=['gas','z_km'],
                    coords={'gas':['O','O2','N2'],
                             'z_km':z},data=glowfort.cglow.sion)


    chdir(str(oldcwd))

    return ver,photIon,isrparam,phitop,zceta,sza,tez,sion





def verprodloss(t,glatlon,flux,EK,makeplot=[None],odir=None,zlim=None):
    """ for a single time, computes VER, production, and loss vs. unit input flux
    inputs:
    -------
    t: a single datetime() when the eigenprofiles should be computed (solar zenith angle computed in Fortran code)
    glatlon: geographic coordinates of site (magnetic coordinates computed in Fortran code)
    flux: a vector of scaled differential number flux for each energy bin, scaled to "unit" energy flux for your eigenprofile input
    E0: a vector of energies [eV] to compute unit responses

    """
    assert isinstance(t,datetime)
    (glat,glon) = glatlon

    vers = None
    for e in EK:
        print('{t} E0: {e:.0f}')

        ver,photIon,isr,phitop,zceta,sza,prate,lrate,tez,sion,UV = runglowaurora(flux,e,t,glat,glon)
     
        if vers is None:
            prates = DataArray(data=np.empty((EK.size,prate.type.size,prate.z_km.size,prate.reaction.size)),
                               dims=['eV','type','z_km','reaction'],
                               coords=[EK,prate.type,prate.z_km,prate.reaction])

            lrates = DataArray(data=np.empty((EK.size,lrate.type.size,lrate.z_km.size,lrate.reaction.size)),
                               dims=['eV','type','z_km','reaction'],
                               coords=[EK,lrate.type, lrate.z_km, lrate.reaction])

            vers = DataArray(data=np.empty((EK.size,ver.z_km.size,ver.wavelength_nm.size)),
                             dims=['eV','z_km','wavelength_nm'],
                             coords=[EK,ver.z_km,ver.wavelength_nm])

            tezs = DataArray(data=np.empty((tez.size,EK.size)),
                             dims=['z_km','eV'],
                             coords=[tez.z_km,EK] )

        prates.loc[e,...] = prate.loc['final',...]
        lrates.loc[e,...] = lrate.loc['final',...]
        vers.loc[e,...]= ver
        tezs.loc[:,e] = tez

        #plotaurora(phitop,ver,flux,sza,zceta,photIon,isr,dtime,glat,glon,e0,zlim,makeplot,odir)


    return vers,photIon,isr,phitop,zceta,sza,prates,lrates,tezs,sion

def ekpcolor(eigen):
    if isinstance(eigen,DataArray):
        e0 = eigen.loc[:,'low'].values
        eEnd = eigen.loc[:,'high'][-1]
        diffnumflux = eigen.loc[:,'flux'].values
    elif isinstance(eigen,(str,Path)):
        eigen = Path(eigen).expanduser()
        if eigen.suffix == '.csv':
            e0 =   np.loadtxt(eigen, usecols=[0], delimiter=',')
            eEnd = np.loadtxt(eigen, usecols=[1], delimiter=',')[-1]
            diffnumflux = None
        elif eigen.suffix == '.h5':
            bins = read_hdf(eigen)
            e0 = bins['low'].values
            eEnd = bins['high'].iloc[-1]
            diffnumflux = bins['flux'].values
        else:
            raise ValueError(f'I do not understand what file you want me to read {eigen}')
    else:
        raise ValueError(f'unknown data type {type(eigen)}')

    return np.append(e0,eEnd),e0,diffnumflux

def makeeigen(EK,diffnumflux,T,glatlon,makeplot=[None],odir=None,zlim=None):
    if isinstance(T,datetime):
        T=[T]

    ver = None
    if len(T) > 1:
       logging.error('more than one time not yet implmented. You will get last time')

    for t in T:
        v,photIon,isr,phitop,zceta,sza,prate,lrate,tez,sion = verprodloss(t,glatlon,diffnumflux,EK, makeplot,odir,zlim)

#TODO time stack
        ver=v; prates=prate; lrates=lrate; tezs=tez

    return ver,photIon,isr,phitop,zceta,sza,prates,lrates,tezs,sion
