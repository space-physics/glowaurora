#!/bin/bash
# Michael Hirsch

#cleanup from old build (force recompilation)
rm -rf build dist glowaurora.egg-info

# call f2py to compile fortran and copy into install directory
python$1 setup.py install

if [[ $? -eq 0 ]]; then
# FIXME manual monkeypatch of fortran data due to missing path
    (
    cd $HOME
    modpath=$(python$1 -c "import glowaurora; print(glowaurora.__path__[0])")

    mv -v $modpath/../glowfort.* $modpath/
    )
else
    echo error in compiling GLOW package
    exit 1
fi
