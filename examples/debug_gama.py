#!/usr/bin/env python3
"""
simulates single float problem in existing Fortran code
"""
from __future__ import division
from numpy import float32,isfinite

Prodwnn1 = float32(-1.15867674E+09)
PROdwn= float32(1.00000000E-30)
Prodwnp1= float32(-640330880.)

PROdup= float32(-998047680.)

Prodn1= float32(0.)
Prod= float32(1.00000000E-30)
Prodp1= float32(0.)

T1= float32(6.20506251E-12)
T2= float32(2.29550547E-08)
Alpha= float32(1.47421190E-06)
Del2= float32(200000.000)


gama = (
        (Prod/float32(2)) * (-T1 - T2 - Alpha - (Prodp1 - Prodn1) /Prod / Del2)
                  + PROdwn
                  * ( -Alpha - T2 - (Prodwnp1 - Prodwnn1) /PROdwn / Del2 )
                  - PROdup * T1
       )

assert isfinite(gama)  #simulates single float problem in existing Fortran code