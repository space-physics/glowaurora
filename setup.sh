#!/bin/bash
# Michael Hirsch

#cleanup from old build (force recompilation)
[[ -z $2 ]] && { echo "removing old fortran .so"; rm -rf build dist glowaurora.egg-info; }

# call f2py to compile fortran and copy into install directory
python$1 setup.py install
# ^^^ must be INSTALL to copy the .dat files and iri file et al to the proper .so relationship

if [[ $? -eq 0 ]] && [[ -z $2 ]]; then
# FIXME manual monkeypatch of fortran data due to missing path
    (
    cd $HOME
    modpath=$(python$1 -c "import glowaurora; print(glowaurora.__path__[0])")
    echo "modpath = $modpath"

    mv -v $modpath/../glowfort.* $modpath/
    )
elif [[ -z $2 ]]; then
    echo '******* error in compiling GLOW package'
    exit 1
fi
