.. image:: https://travis-ci.org/scienceopen/glowaurora.svg
    :target: https://travis-ci.org/scienceopen/glowaurora

.. image:: https://coveralls.io/repos/scienceopen/glowaurora/badge.svg?branch=master&service=github 
    :target: https://coveralls.io/github/scienceopen/glowaurora?branch=master 
    
.. image:: https://codeclimate.com/github/scienceopen/glowaurora/badges/gpa.svg
   :target: https://codeclimate.com/github/scienceopen/glowaurora
   :alt: Code Climate

=============
glow-aurora
=============
`Stan Solomon's  GLOW Auroral model <http://download.hao.ucar.edu/pub/stans/glow/>`_ -- now in Python!

:Fortran author: Stan Solomon
:Python API author: Michael Hirsch

.. contents::

.. image:: examples/demo_out.png
   :alt: vertical profiles of VER

.. image:: examples/demo_in.png
   :alt: diff num flux input

Installation
============
Note, if you don't already have Numpy installed, just run the script below twice. 
First time installs Numpy if you don't have it, and then the main program.::

   python setup.py develop

Examples
========

Self-test f2py
--------------
This self-test should give zero errors. This tests the Fortran code from Python.::
  
  ./test/test.py -v


volume emission rate plots 
--------------------------
To produce the plots seen at the Github site::

  python examples/demo_glowaurora.py

with the volume emission rate and intermediate
processes modeled for the given primary electron precipitation input. You can make
this more generally useful as eigenprofiles in the next section.

production/loss rate eigenprofiles
----------------------------------
This requires two steps:

1. Generate unit input differential number flux vs. energy
2. Compute ionospheric energy deposition and hence production/loss rates for the modeled kinetic chemistries (12 in total)

This is handled by the script ``gridaurora/MakeIonoEigenprofile.py``

Papers
======
(Thanks to Stephen Kaeppler to pointing these out)

http://download.hao.ucar.edu/pub/stans/papers/BaileyJGR2002.pdf

http://download.hao.ucar.edu/pub/stans/papers/SolomonJGR1988.pdf

Appendix (Not necessary for the typical user)
=============================================
has been tested at various times with gfortran 4.8 - 5.2 on Windows and Linux.
Any issues please contact me. I will try to make my best
effort, several researchers are using this code already.


Download the GLOW v0.973 source code from Stan Solomon
------------------------------------------------------
Stan's team has been working on a new version, Modern Fortran, looked beautiful
from a sneak peek, but for now we'll be satiated with the original.::

  wget -r -np -nc -nH --cut-dirs=4 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/glow/v0.973/

Download Stan's copy of IRI files
---------------------------------
Stan tweaked IRI90 slightly, here's the copy he uses.::

  wget -r -np -nc -nH --cut-dirs=3 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/iri/


compile the Fortran code by itself
----------------------------------
The Fortran program used by itself spits out a lot of text as its output::

  cd fortran
  make

F2PY compile the Fortran code for use from Python
-------------------------------------------------
::

   f2py3 -m glowfort -c egrid.f maxt.f glow.f vquart.f gchem.f ephoto.f solzen.f rcolum.f etrans.f exsect.f ssflux.f snoem.f snoemint.f geomag.f nrlmsise00.f qback.f fieldm.f iri90.f aurora_sub.f --quiet

You can pick a specific compiler by adding the ``--f90exec=`` option. For example
you could use the Intel Fortran Compilerseparately for this by starting with::

    f2py3 --f90exec=gfortran-5.2 -m glowfort

and so on.


Fortran self-test
-----------------
::

  ./auroraexample < aurexample.in > aurtest.out

observe that aurtest.out is almost exactly equal to reference/aurexample.out, to the least digit of precision.


Licensing
=========
original Fortran code in directory ``fortran/`` as obtained from http://download.hao.ucar.edu/pub/stans/glow/:

"This software is part of the GLOW model.  Use is governed by the Open Source Academic Research License
Agreement contained in the file glowlicense.txt."


Python code and modifications to original Fortran code:  GNU Affero GPLv3+
