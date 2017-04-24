#!/usr/bin/env python
"""
Example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
"""
from __future__ import print_function
try:
    from pathlib import Path
    Path().expanduser()
except (ImportError,AttributeError):
    from pathlib2 import Path
#
import logging
from sys import stderr
from datetime import datetime
from xarray import DataArray
from numpy import hstack,degrees,zeros_like,ndarray,atleast_1d,empty,append,loadtxt
from pandas import read_hdf
from os import chdir
#
from sciencedates import datetime2yd, find_nearest
from gridaurora import readmonthlyApF107
from gridaurora.zglow import glowalt
import glowfort
glowpath=Path(__file__).parents[1]
oldcwd = Path.cwd()

def runglowaurora(eflux,e0,t0,glat,glon,f107apfn=None,f107a=None,f107=None,f107p=None,ap=None):
#%% (-2) check/process user inputs
    assert isinstance(eflux,(float,int,ndarray))
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

    eflux = atleast_1d(eflux)

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

    lamb=[3371., 4278., 5200., 5577., 6300., 7320., 10400., 3466., 7774., 8446., 3726.]
    products=['nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+','nO','nO2','nN2','nNO']

    ver = DataArray(dims=['z_km','wavelength_nm'],
                    coords={'z_km':z,
                            'wavelength_nm':lamb},
                    data=zeta[:,:11])
    photIon = DataArray(dims=['z_km','type'],
                        coords={'z_km':z,
                                'type':['photoIoniz','eImpactIoniz','ne']+products},
                   data=hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)))

    isrparam = DataArray(dims=['z_km','param'],
                        coords={'z_km':z,'param':['ne','Te','Ti']},
                         data=isr)

    phitop = DataArray(coords={'eV':phi[:,0]},
                       data=phi[:,2])

    zceta = DataArray(dims=['z_km','wavelength_nm','type'],
                    data=glowfort.cglow.zceta.T[:,:11,:])  #Nalt x Nwavelengths  xNproductionEmissions
    zceta['z_km']=z
    zceta['wavelength_nm']=lamb
    #FIxmE what is 3rd axis?


    sza = degrees(glowfort.cglow.sza)

    tez = DataArray(coords={'z_km':z},
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
        print('{} E0: {:.0f}'.format(t,e))

        ver,photIon,isr,phitop,zceta,sza,prate,lrate,tez,sion = runglowaurora(flux,e,t,glat,glon)
        if vers is None:
            prates = DataArray(data=empty((EK.size,prate.type.size,prate.z_km.size,prate.reaction.size)),
                               dims=['eV','type','z_km','reaction'],
                               coords=[EK,prate.type,prate.z_km,prate.reaction])

            lrates = DataArray(data=empty((EK.size,lrate.type.size,lrate.z_km.size,lrate.reaction.size)),
                               dims=['eV','type','z_km','reaction'],
                               coords=[EK,lrate.type, lrate.z_km, lrate.reaction])

            vers = DataArray(data=empty((EK.size,ver.z_km.size,ver.wavelength_nm.size)),
                             dims=['eV','z_km','wavelength_nm'],
                             coords=[EK,ver.z_km,ver.wavelength_nm])

            tezs = DataArray(data=empty((tez.size,EK.size)),
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
            e0 =   loadtxt(str(eigen),usecols=[0],delimiter=',')
            eEnd = loadtxt(str(eigen),usecols=[1],delimiter=',')[-1]
            diffnumflux = None
        elif eigen.suffix == '.h5':
            bins = read_hdf(str(eigen))
            e0 = bins['low'].values
            eEnd = bins['high'].iloc[-1]
            diffnumflux = bins['flux'].values
        else:
            raise ValueError('I do not understand what file you want me to read {}'.format(eigen))
    else:
        raise ValueError('unknown data type {}'.format(type(eigen)))

    return append(e0,eEnd),e0,diffnumflux

def makeeigen(EK,diffnumflux,T,glatlon,makeplot=[None],odir=None,zlim=None):
    if isinstance(T,datetime):
        T=[T]

    ver = None
    if len(T) > 1:
        print('more than one time not yet implmented. You will get last time',file=stderr)

    for t in T:
        v,photIon,isr,phitop,zceta,sza,prate,lrate,tez,sion = verprodloss(t,glatlon,diffnumflux,EK, makeplot,odir,zlim)

#TODO time stack
        ver=v; prates=prate; lrates=lrate; tezs=tez

    return ver,photIon,isr,phitop,zceta,sza,prates,lrates,tezs,sion
