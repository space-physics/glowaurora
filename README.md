[![image](https://zenodo.org/badge/34395725.svg)](https://zenodo.org/badge/latestdoi/34395725)
[![image](https://travis-ci.org/scivision/glowaurora.svg)](https://travis-ci.org/scivision/glowaurora)
[![image](https://coveralls.io/repos/github/scivision/glowaurora/badge.svg?branch=master)](https://coveralls.io/github/scivision/glowaurora?branch=master)
[![Maintainability](https://api.codeclimate.com/v1/badges/0345889286d7ae6a8307/maintainability)](https://codeclimate.com/github/scivision/glowaurora/maintainability)

New [Modern Fortran GLOW](https://github.com/scivision/GLOW) with Python, Matlab, or MPI.

# GLOW 0.973

GLOW 0.973 Airglow/Auroral [model](http://download.hao.ucar.edu/pub/stans/glow/) -- now in Python!

![incident energy and VER](examples/ver.png)

![vertical profiles](examples/demo_out.png)

![diff num flux input](examples/demo_in.png)

## Install
A Fortran compiler is required:

* Linux: `apt install gfortran`
* MacOS: `brew install gcc`
* [Windows](https://www.scivision.co/windows-gcc-gfortran-cmake-make-install)
 


```sh
pip install -e .
```

## Examples

### Self-test f2py

This self-test should give zero errors. 
This tests the Fortran code from Python.:

    pytest -v

### volume emission rate plots

To produce the plots seen at the Github site:

    python RunGlow.py

with options including:

* `-t --simtime` time of simulation 
* `-c --latlon`  geographic coordinate (lat,lon) [degrees] 
* `-q --flux`  total flux

with the volume emission rate and intermediate processes modeled for the
given primary electron precipitation input. 
You can make this more
generally useful as eigenprofiles in the next section.

### production/loss rate eigenprofiles

This requires two steps:

1.  Generate unit input differential number flux vs. energy
2.  Compute ionospheric energy deposition and hence production/loss
    rates for the modeled kinetic chemistries (12 in total)

This is handled by the script `gridaurora/MakeIonoEigenprofile.py`

### Matlab access to Glow

Matlab can use Glow via the Python interface, as in the example
`glow.m`.

## Papers

(Thanks to Stephen Kaeppler to pointing these out)

* [Bailey 2002](http://download.hao.ucar.edu/pub/stans/papers/BaileyJGR2002.pdf)
* [Solomon 1988](http://download.hao.ucar.edu/pub/stans/papers/SolomonJGR1988.pdf)

## Appendix (Not necessary for the typical user)

### Download the GLOW v0.973 source code from Stan Solomon

Stan's team has released GLOW v0.98 using modern Fortran, but here's
the original.

```sh
wget -r -np -nc -nH --cut-dirs=4 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/glow/v0.973/
```

### Download Stan's copy of IRI files

Stan tweaked IRI90 slightly, here's the copy he uses.

```sh
wget -r -np -nc -nH --cut-dirs=3 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/iri/
```

### compile the Fortran code by itself

The Fortran program used by itself spits out a lot of text as its
output:

```sh
cd bin
cmake ../src
cmake --build .
```

### Fortran self-test

Auroral example

```sh
./auroraexample < aurexample.in > aurtest.dat
```

High energy example

```sh
./hexexample < hexexample.in > hextest.dat
```

## Notes


### Licensing

original Fortran code in directory `fortran/` as obtained from
<http://download.hao.ucar.edu/pub/stans/glow/>: "This software is part
of the GLOW model. Use is governed by the Open Source Academic Research
License Agreement contained in the file glowlicense.txt."

Python code and modifications to original Fortran code: GNU Affero GPLv3+
