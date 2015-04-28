# glow-aurora 
[Stan Solomon's GLOW Auroral model](http://download.hao.ucar.edu/pub/stans/glow/) -- now in Python!

![Aurora VER demo](http://blogs.bu.edu/mhirsch/files/2015/04/plotglow_panel.png)
 
Installation:
-------------
```
f2py --opt='-fno-align-commons' -m aurora -c glow.f gchem.f ephoto.f solzen.f rcolum.f etrans.f exsect.f vquart.f ssflux.f snoem.f90 snoemint.f geomag.f nrlmsise00.f qback.f fieldm.f iri90.f aurora_sub.f
```
if you want to use the GLOW gridder, separately and additionally do:
```
f2py -m glowgrid -c egrid.f90 maxt.f90
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
