#!/usr/bin/env python3
"""
Trivial example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
bostonmicrowave.com
"""
from matplotlib.pyplot import figure, show
from pandas import DataFrame
from datetime import datetime
from dateutil.parser import parse
#
from fortrandates import datetime2gtd
try:
    from aurora import aurora
except ImportError as e:
    exit('you must compile with f2py first. See README.md')

def demoaurora(iyd,utsec):
    z,zeta = aurora(iyd,utsec,70,0,100,100,100,4,1,1000)
    ver = DataFrame(index=z, data=zeta.T[:,:10],
                    columns=[3371, 4278, 5200, 5577, 6300,7320,10400,3466,
                             7774, 8446])
    return ver

def plotaurora(ver,dtime):
    ax = figure().gca()
    ax.plot(ver.values,ver.index)
    ax.set_xlabel('VER')
    ax.set_ylabel('altitude [km]')
    ax.grid(True)
    ax.legend(ver.columns)
    ax.set_title(str(dtime))

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="Stan Solomon's GLOW auroral model")
    p.add_argument('simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',type=str,nargs='?',default=None)
    p = p.parse_args()

    if p.simtime is None:
        dtime = datetime.now()
    else:
        dtime = parse(p.simtime)


    iyd,utsec = datetime2gtd(dtime)[:2]

    ver = demoaurora(iyd,utsec)
    plotaurora(ver,dtime)
    show()