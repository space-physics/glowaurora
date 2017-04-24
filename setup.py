#!/usr/bin/env python
req = ['nose','python-dateutil','numpy','pandas','xarray','matplotlib','seaborn','astropy','h5py','cython','pathlib2']
pipreq=['pymap3d','sciencedates','gridaurora']
# %%
import pip
try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception:
    pip.main(['install'] + req)
pip.main(['install'] + pipreq)
# %%
import setuptools #needed to enable develop
from numpy.distutils.core import setup,Extension
from glob import glob
from os.path import join

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


# %% 
setup(name='glowaurora',
      packages=['glowaurora'],
      version='1.0.2',
      author='Michael Hirsch, Ph.D.',
      description='Model of auroral and airglow emissions',
      url='https://github.com/scivision/glowaurora',
#      package_dir={'glowaurora': 'glowaurora'}, #not working
 #     package_data={'glowaurora': ['fortran/*.dat']}, #not working, use data_files
      include_package_data=True,
      ext_modules=ext,
      data_files=[('glowaurora',fortdata),
                  ('glowaurora/iri',iridata)
                  ], #must have data_files to copy *.dat to site-packages
      )


