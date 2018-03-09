#!/usr/bin/env python
install_requires = ['python-dateutil','numpy','xarray',
                   'sciencedates','gridaurora']
tests_require=['pytest','nose','coveralls']
# %%
from setuptools import find_packages
from numpy.distutils.core import setup,Extension
from glob import glob
from os.path import join

# f2py -m aurora -c egrid.f maxt.f glow.f vquart.f gchem.f ephoto.f solzen.f rcolum.f etrans.f exsect.f ssflux.f snoem.f snoemint.f geomag.f nrlmsise00.f qback.f fieldm.f aurora_sub.f

fortranfiles=['egrid.f','maxt.f','glow.f','rout.f',
              'vquart.f','gchem.f','ephoto.f','solzen.f','rcolum.f',
              'etrans.f','exsect.f','ssflux.f','snoem.f','snoemint.f',
              'geomag.f','qback.f','fieldm.f',
              'nrlmsise00.f','iri90.f',
              'aurora_sub.f','dayglow_sub.f']

root='fortran'

fortranpaths = [join(root,f) for f in fortranfiles]
fortdata = glob(join(root,'*.dat'))
iridata = glob(join('iri','*.asc'))
#%% prelim
ext=[Extension(name='glowfort',
               sources=fortranpaths,
               f2py_options=['--quiet'],
               #extra_f77_compile_args=['-O0'],
               #['-finit-local-zero'] #not needed
)]
               #include_dirs=[root],
               #library_dirs=[root])]


# %%
setup(name='glowaurora',
      packages=find_packages(),
      version='1.2.1',
      author='Michael Hirsch, Ph.D.',
      description='Model of auroral and airglow emissions',
      long_description=open('README.rst').read(),
      url='https://github.com/scivision/glowaurora',
#      package_dir={'glowaurora': 'glowaurora'}, #not working
 #     package_data={'glowaurora': ['fortran/*.dat']}, #not working, use data_files
      include_package_data=True,
      ext_modules=ext,
      data_files=[('glowaurora',fortdata),
                  ('glowaurora/iri',iridata)
                  ], #must have data_files to copy *.dat to site-packages
      install_requires=install_requires,
      python_requires='>=3.6',
      tests_require=tests_require,
      extras_require={'plot':['matplotlib>=2.2','seaborn'],
                        'io':['h5py','astropy','pandas',],
                      'tests':tests_require},
        classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
       'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3.6',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      ],
      )

