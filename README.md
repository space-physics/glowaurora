[![Code Climate](https://codeclimate.com/github/scienceopen/glowaurora/badges/gpa.svg)](https://codeclimate.com/github/scienceopen/glowaurora)

# glow-aurora 
[Stan Solomon's  GLOW Auroral model](http://download.hao.ucar.edu/pub/stans/glow/) -- now in Python!

![Aurora VER demo](demo_out.png)
![flux input](demo_in.png)

Caution:
--------
It is currently (0.973) necessary to restart the Python kernel for each simulation run. This is due to the "save" statement in ssflux.f.
With the refactoring of the code underway to Fortran90, it is hoped this blanket save statement can be eliminated.

For safety's sake, run this program from the command line (Terminal) to ensure you get a fresh import (flushing all variables).
 
Installation:
-------------
```
git clone --recursive https://github.com/scienceopen/glow-aurora

make -f Makefile.f2py
```


Note: If on Windows and using MinGW compiler, add the option ``` --compiler=mingw32 ```

Yes, even though [you're using 64-bit MinGW](http://blogs.bu.edu/mhirsch/2015/04/f2py-running-fortran-code-in-python-on-windows/)

```
python demo_aurora.py
```
will show modeled VER vs. altitude for the input parameter set.


Papers:
------
(Thanks to Stephen Kaeppler to pointing these out)

http://download.hao.ucar.edu/pub/stans/papers/BaileyJGR2002.pdf 

http://download.hao.ucar.edu/pub/stans/papers/SolomonJGR1988.pdf

Appendix (Not necessary for the typical user):
----------------------------------------------
### Download the latest source code from Stan Solomon:
``` 
wget -r -np -nc -nH --cut-dirs=4 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/glow/v0.973/
```

### Download Stan's copy of IRI files:
```
wget -r -np -nc -nH --cut-dirs=3 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/iri/
```

### compile the Fortran code by itself
```
make
```

### Fortran self-test
```
./auroraexample < aurexample.in > aurtest.out
```
observe that aurtest.out is almost exactly equal to aurexample.out, to the least digit of precision.
