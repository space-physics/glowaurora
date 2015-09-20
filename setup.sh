#!/bin/bash
python$1 setup.py install

# FIXME manual monkeypatch of fortran data due to missing path
(
cd $HOME
modpath=$(python$1 -c "import glowaurora; print(glowaurora.__path__[0])")

mv -v $modpath/../glowfort.* $modpath/
)

