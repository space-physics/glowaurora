#!/usr/bin/env python3
"""
loads and plots eigenprofiles vs time

"""
from __future__ import division,absolute_import
import h5py
from os.path import expanduser
from datetime import datetime
from pytz import UTC
#
from histfeas.plotsnew import ploteig
from transcarread.readTranscar import SimpleSim
from RunLoop import ekpcolor

from argparse import ArgumentParser
p = ArgumentParser()
p.add_argument('-E','--eigenfn',help='generate eigenprofiles using energies in this csv file',default='~/code/transcar/transcar/BT_E1E2prev.csv')
p.add_argument('-i','--infn',help='input hdf5 to read',default='~/data/transcareigen.h5')
p = p.parse_args()


EKpcolor,e0 = ekpcolor(p.eigenfn)

sim = SimpleSim(filt='none',inpath=None,reacreq='')

with h5py.File(expanduser(p.infn),'r') as f:
    if len(f['/eigenprofile'].shape) == 3:
        for i,v in enumerate(f['/eigenprofile']): #for each time
            try:
                t = datetime.fromtimestamp(f['ut1_unix'].value[i],tz=UTC)
            except:
                t=None
            ploteig(EKpcolor,f['/altitude'].value,v,(None,)*6,sim,t)
    elif len(f['/eigenprofile'].shape) ==2: #old single time
        ploteig(EKpcolor,f['/altitude'].value,f['/eigenprofile'].value,(None,)*6,sim)