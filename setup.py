#!/usr/bin/env python

import setuptools  # noqa: F401
from numpy.distutils.core import setup, Extension
from glob import glob
from os.path import join
from pathlib import Path
import os


if os.name == 'nt':
    sfn = Path(__file__).parent / 'setup.cfg'
    stxt = sfn.read_text()
    if '[build_ext]' not in stxt:
        with sfn.open('a') as f:
            f.write("[build_ext]\ncompiler = mingw32")

# f2py -m aurora -c egrid.f maxt.f glow.f vquart.f gchem.f ephoto.f solzen.f
# rcolum.f etrans.f exsect.f ssflux.f snoem.f snoemint.f geomag.f nrlmsise00.f qback.f fieldm.f aurora_sub.f

fortranfiles = ['egrid.f', 'maxt.f', 'glow.f', 'rout.f',
                'vquart.f', 'gchem.f', 'ephoto.f', 'solzen.f', 'rcolum.f',
                'etrans.f', 'exsect.f', 'ssflux.f', 'snoem.f', 'snoemint.f',
                'geomag.f', 'qback.f', 'fieldm.f',
                'nrlmsise00.f', 'iri90.f',
                'aurora_sub.f', 'dayglow_sub.f']

root = 'src'

fortranpaths = [join(root, f) for f in fortranfiles]
fortdata = glob(join(root, '*.dat'))
iridata = glob(join('iri', '*.asc'))
# %% prelim
ext = [Extension(name='glowfort',
                 sources=fortranpaths,
                 f2py_options=[],
                 # extra_f77_compile_args=['-O0'],
                 # ['-finit-local-zero'] #not needed
                 )]

setup(  # package_dir={'glowaurora': 'glowaurora'}, #not working
    #     package_data={'glowaurora': ['fortran/*.dat']}, #not working, use data_files
    ext_modules=ext,
    data_files=[('glowaurora', fortdata),
                ('glowaurora/iri', iridata)
                ],  # must have data_files to copy *.dat to site-packages
)
