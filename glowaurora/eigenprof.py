from pathlib import Path
from warnings import warn
from datetime import datetime
from numpy import loadtxt,append,empty
from pandas import read_hdf
from xarray import DataArray
#
from .runglow import runglowaurora

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
        warn('more than one time not yet implmented. You will get last time')

    for t in T:
        v,photIon,isr,phitop,zceta,sza,prate,lrate,tez,sion = verprodloss(t,glatlon,diffnumflux,EK, makeplot,odir,zlim)

#TODO time stack
        ver=v; prates=prate; lrates=lrate; tezs=tez

    return ver,photIon,isr,phitop,zceta,sza,prates,lrates,tezs,sion
