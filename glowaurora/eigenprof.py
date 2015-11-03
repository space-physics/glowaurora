from __future__ import division,absolute_import
from numpy import loadtxt,append
from pandas import DataFrame,Panel,Panel4D,read_hdf
from os.path import expanduser
#
from .runglow import runglowaurora

def verprodloss(t,glatlon,flux,EK,makeplot,odir,zlim):
    """ for a single time, computes VER, production, and loss vs. unit input flux
    inputs:
    -------
    t: a single datetime() when the eigenprofiles should be computed (solar zenith angle computed in Fortran code)
    glatlon: geographic coordinates of site (magnetic coordinates computed in Fortran code)
    flux: a vector of scaled differential number flux for each energy bin, scaled to "unit" energy flux for your eigenprofile input
    E0: a vector of energies [eV] to compute unit responses

    """

    (glat,glon) = glatlon

    vers = None #sentinal
    for e in EK:
        print('{} E0: {:.0f}'.format(t,e))

        ver,photIon,isr,phitop,zceta,sza,prate,lrate,tez,sion = runglowaurora(flux,e,t,glat,glon)
        if vers is None:
            prates=Panel(items=EK, major_axis=prate.major_axis, minor_axis=prate.minor_axis)
            lrates=Panel(items=EK, major_axis=lrate.major_axis, minor_axis=lrate.minor_axis)
            vers = Panel(items=EK, major_axis=ver.index,        minor_axis=ver.columns)
            tezs = DataFrame(index=ver.index,columns=EK)

        vers[e] = ver
        prates[e]=prate['final']
        lrates[e]=lrate['final']
        tezs[e]=tez

        #plotaurora(phitop,ver,flux,sza,zceta,photIon,isr,dtime,glat,glon,e0,zlim,makeplot,odir)


    return vers,photIon,isr,phitop,zceta,sza,prates,lrates,tezs,sion

def ekpcolor(eigen):
    if isinstance(eigen,DataFrame):
        e0 = eigen['low'].values
        eEnd = eigen['high'].iloc[-1]
        diffnumflux = eigen['flux'].values
    else:
        if eigen.endswith('.csv'):
            e0 =   loadtxt(expanduser(eigen),usecols=[0],delimiter=',')
            eEnd = loadtxt(expanduser(eigen),usecols=[1],delimiter=',')[-1]
            diffnumflux = None
        elif eigen.endswith('.h5'):
            bins = read_hdf(expanduser(eigen))
            e0 = bins['low'].values
            eEnd = bins['high'].iloc[-1]
            diffnumflux = bins['flux'].values
        else:
            raise ValueError('I do not understand what file you want me to read {}'.format(eigen))

    return append(e0,eEnd),e0,diffnumflux

def makeeigen(EK,diffnumflux,T,glatlon,makeplot,odir,zlim):
    ver = None

    for t in T:
        v,photIon,isr,phitop,zceta,sza,prate,lrate,tez,sion = verprodloss(t,glatlon,diffnumflux,EK, makeplot,odir,zlim)
        if ver is None:
            ver =  Panel4D(labels=T,items=v.items,    major_axis=v.major_axis,    minor_axis=v.minor_axis)
            prates=Panel4D(labels=T,items=prate.items,major_axis=prate.major_axis,minor_axis=prate.minor_axis)
            lrates=Panel4D(labels=T,items=lrate.items,major_axis=lrate.major_axis,minor_axis=lrate.minor_axis)
            tezs=Panel(items=T,major_axis=tez.index,minor_axis=tez.columns)

        ver[t] = v # v is a 3-D Panel
        prates[t]=prate
        lrates[t]=lrate
        tezs[t]=tez

    return ver,photIon,isr,phitop,zceta,sza,prates,lrates,tezs,sion