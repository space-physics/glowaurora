#!/usr/bin/env python
import setuptools #needed to enable develop
from glob import glob
from os.path import join
from numpy.distutils.core import setup,Extension

# f2py -m aurora -c egrid.f maxt.f glow.f vquart.f gchem.f ephoto.f solzen.f rcolum.f etrans.f exsect.f ssflux.f snoem.f snoemint.f geomag.f nrlmsise00.f qback.f fieldm.f aurora_sub.f

# FUTURE: iri90.f not included in f2py b/c we use the pyiri90 package
# FUTURE: nrlmsise00.f not included .... msise00 package
fortranfiles=['egrid.f','maxt.f','glow.f',
              'vquart.f','gchem.f','ephoto.f','solzen.f','rcolum.f',
              'etrans.f','exsect.f','ssflux.f','snoem.f','snoemint.f',
              'geomag.f','qback.f','fieldm.f',
              'nrlmsise00.f','iri90.f',
              'aurora_sub.f']

root='fortran'

fortranpaths = [join(root,f) for f in fortranfiles]
fortdata = glob(join(root,'*.dat'))
iridata = glob(join('iri','*.asc')) #in pyiri90
#%% prelim
ext=[Extension(name='glowfort',
               sources=fortranpaths,
               f2py_options=['--quiet'],
               #extra_f77_compile_args=['-finit-local-zero'] #not needed
)]
               #include_dirs=[root],
               #library_dirs=[root])]

req = ['nose','python-dateutil','numpy','pandas','xarray','matplotlib','seaborn','astropy','h5py','cython',
        'pymap3d','sciencedates','gridaurora']

#%% install
setup(name='glowaurora',
      packages=['glowaurora'],
      author='Michael Hirsch, Ph.D.',
#      package_dir={'glowaurora': 'glowaurora'}, #not working
 #     package_data={'glowaurora': ['fortran/*.dat']}, #not working, use data_files
      include_package_data=True,
      ext_modules=ext,
      data_files=[('glowaurora',fortdata),
                  ('glowaurora/iri',iridata)
                  ], #must have data_files to copy *.dat to site-packages

	  install_requires=req,
      dependency_links = [
          'https://github.com/scivision/gridaurora/tarball/master#egg=gridaurora',
                            ],
      )


