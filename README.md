# glow-aurora
Stan Solomon's GLOW Auroral model -- now in Python!

Installation:
-------------
```
f2py3 -m aurora -c aurora_sub.f maxt.f glow.f gchem.f ephoto.f solzen.f rcolum.f etrans.f exsect.f vquart.f ssflux.f egrid.f snoemint.f snoem.f geomag.f nrlmsise00.f qback.f fieldm.f iri90.f
```

```
python
import aurora
```


Appendix (Not necessary for the typical user):
----------------------------------------------
### Download the latest source code from Stan Solomon:
``` 
wget -r -np -nc -nH --cut-dirs=4 --random-wait --wait 1 -R "index.html*" http://download.hao.ucar.edu/pub/stans/glow/v0.973/
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
