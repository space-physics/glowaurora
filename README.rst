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
::

   python setup.py develop

Examples
========

Self-test f2py
--------------
::
  
  cd test
  nosetests -v test.py

the self-test should give zero errors

volume emission rate plots 
--------------------------
::

  python examples/demo_glowaurora.py

this produces the plots seen here, with the volume emission rate and intermediate
processes modeled for the given primary electron precipitation input.

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
has been tested at various times with gfortran 4.6 - 5.2 on Windows and Linux.
Any issues please contact me.


Download the GLOW v0.973 source code from Stan Solomon
-------------------------------------------------
::

  wget -r -np -nc -nH --cut-dirs=4 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/glow/v0.973/

Download Stan's copy of IRI files
---------------------------------
::

  wget -r -np -nc -nH --cut-dirs=3 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/iri/


compile the Fortran code by itself
----------------------------------
::

  cd fortran
  make

F2PY compile the Fortran code for use from Python
-------------------------------------------------
::

 f2py -m glowfort -c egrid.f maxt.f glow.f vquart.f gchem.f ephoto.f solzen.f rcolum.f etrans.f exsect.f ssflux.f snoem.f snoemint.f geomag.f nrlmsise00.f qback.f fieldm.f iri90.f aurora_sub.f --quiet


Fortran self-test
-----------------
::

  ./auroraexample < aurexample.in > aurtest.out

observe that aurtest.out is almost exactly equal to reference/aurexample.out, to the least digit of precision.

Windows operating system
========================
On Windows, consider `factors like <https://scivision.co/f2py-running-fortran-code-in-python-on-windows/>`_


Licensing
=========
original Fortran code in directory ``fortran/`` as obtained from http://download.hao.ucar.edu/pub/stans/glow/:

"This software is part of the GLOW model.  Use is governed by the Open Source Academic Research License
Agreement contained in the file glowlicense.txt."


Python code and modifications to original Fortran code:  GNU Affero GPLv3+
