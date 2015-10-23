from __future__ import division,absolute_import
from numpy import loadtxt,append
from pandas import DataFrame,Panel,read_hdf
from os.path import expanduser
#
from .runglow import runglowaurora

def verprodloss(t,glatlon,flux,EK,f107a,f107,f107p,ap,makeplot,odir,zlim):
    """ for a single time, computes VER, production, and loss vs. unit input flux
    inputs:
    -------
    t: a single datetime() when the eigenprofiles should be computed (solar zenith angle computed in Fortran code)
    glatlon: geographic coordinates of site (magnetic coordinates computed in Fortran code)
    flux: a vector of scaled differential number flux for each energy bin, scaled to "unit" energy flux for your eigenprofile input
    E0: a vector of energies [eV] to compute unit responses

    """

    (glat,glon) = glatlon

    DFver = DataFrame(); prates=None; lrates=None
    for e in EK:
        print('{} E0: {:.0f}'.format(t,e))

        ver,photIon,isr,phitop,zceta,sza,prate,lrate = runglowaurora(flux,e,
                                                                  t,glat,glon,
                                                                  f107a,f107,f107p,ap)
        if prates is None:
            prates=Panel(items=EK, major_axis=prate.major_axis, minor_axis=prate.minor_axis)
            lrates=Panel(items=EK, major_axis=lrate.major_axis, minor_axis=lrate.minor_axis)

        prates[e]=prate['final']
        #plotaurora(phitop,ver,flux,sza,zceta,photIon,isr,dtime,glat,glon,e0,zlim,makeplot,odir)

        DFver[e] = ver.sum(axis=1)

    return DFver,photIon,isr,phitop,zceta,sza,prates,lrates

def ekpcolor(eigenfn):
    if eigenfn.endswith('.csv'):
        e0 =   loadtxt(expanduser(eigenfn),usecols=[0],delimiter=',')
        eEnd = loadtxt(expanduser(eigenfn),usecols=[1],delimiter=',')[-1]
        diffnumflux = None
    elif eigenfn.endswith('.h5'):
        bins = read_hdf(expanduser(eigenfn))
        e0 = bins['low']
        eEnd = bins['high'].iloc[-1]
        diffnumflux = bins['flux']
    else:
        raise ValueError('I do not understand what file you want me to read {}'.format(eigenfn))

    return append(e0,eEnd),e0,diffnumflux

def makeeigen(eigenfn,T,glatlon,f107a,f107,f107p,ap,makeplot,odir,zlim):
    makeplot.append('eig')
    EKpcolor,EK,diffnumflux = ekpcolor(eigenfn)

    ver = None

    for t in T:
        v,photIon,isr,phitop,zceta,sza,prates,lrates = verprodloss(t,glatlon,diffnumflux,EK,
                                                                   f107a,f107,f107p,ap,
                                                                   makeplot,odir,zlim)
        if ver is None:
            ver = Panel(items=T,major_axis=v.index,minor_axis=v.columns)

        ver.loc[t,:,:] = v

    return ver,photIon,isr,phitop,zceta,sza,EKpcolor,prates,lrates