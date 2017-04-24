.. image:: https://zenodo.org/badge/34395725.svg
   :target: https://zenodo.org/badge/latestdoi/34395725
.. image:: https://travis-ci.org/scivision/glowaurora.svg
    :target: https://travis-ci.org/scivision/glowaurora
.. image:: https://coveralls.io/repos/github/scivision/glowaurora/badge.svg?branch=master
:target: https://coveralls.io/github/scivision/glowaurora?branch=master

    

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

Prereq
======

Linux / BSD / Windows::

    apt install gfortran

Mac::
    
    brew install gcc

Installation
============
::

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

  ./RunGlow.py

with options including:

-t, --simtime   time of simulation (ending in Z for UTC)
-c, --latlon    geographic coordinate (lat,lon) [degrees]
-q, --flux      total flux

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

  cd bin
  cmake ../fortran
  make


Fortran self-test
-----------------
Auroral example::

  ./auroraexample < aurexample.in > aurtest.dat


High energy example::

  ./hexexample < hexexample.in > hextest.dat



Notes
=====

Windows
-------
Strongly suggest using Linux Subsystem for Windows, if you aren't using it yet you should be for your development work.


Licensing
=========
original Fortran code in directory ``fortran/`` as obtained from http://download.hao.ucar.edu/pub/stans/glow/:

"This software is part of the GLOW model.  Use is governed by the Open Source Academic Research License
Agreement contained in the file glowlicense.txt."


Python code and modifications to original Fortran code:  GNU Affero GPLv3+
