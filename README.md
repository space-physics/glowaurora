# glow-aurora
Stan Solomon's Aurora model -- now in Python!




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

