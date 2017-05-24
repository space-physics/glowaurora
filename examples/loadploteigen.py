#!/usr/bin/env python
"""
loads and plots eigenprofiles vs time

"""
from pathlib import Path
import h5py
from datetime import datetime
from pytz import UTC
#
from histfeas.plotsnew import ploteig  # obtained from https://github.com/scivision/histfeas
from transcarread import SimpleSim  # obtained from https://github.com/scivision/transcarread
from glowaurora.eigenprof import ekpcolor

#%% main program
from argparse import ArgumentParser
p = ArgumentParser()
p.add_argument('-E','--eigenfn',help='generate eigenprofiles using energies in this csv file',default='../../gridaurora/zettflux.csv')
p.add_argument('-i','--infn',help='input hdf5 to read',default='~/data/transcareigen.h5')
p = p.parse_args()

fn = Path(p.infn).expanduser()

EKpcolor,e0,diffnumflux = ekpcolor(p.eigenfn)

sim = SimpleSim(filt='none',inpath=None,reacreq='')
#%% load and plot
with h5py.File(str(fn),'r',libver='latest') as f:
    if len(f['/eigenprofile'].shape) == 3:
        for i,v in enumerate(f['/eigenprofile']): #for each time
            try:
                t = datetime.fromtimestamp(f['ut1_unix'][i],tz=UTC)
            except IndexError:
                t=None
            ploteig(EKpcolor,f['/altitude'].value,v,(None,)*6,sim,t)
    elif len(f['/eigenprofile'].shape) ==2: #old single time
        ploteig(EKpcolor,f['/altitude'].value,f['/eigenprofile'].value,(None,)*6,sim)
