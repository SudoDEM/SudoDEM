#!/bin/bash

OS="ubuntu-20.04"

cd ..
mkdir -p 3rdlib/3rdlibs/py

# download pre-compiled tools
wget https://github.com/SudoDEM/3rdlibs/releases/latest/download/tools-${OS}.tar.xz
tar xf tools-${OS}.tar.xz
# download pre-compiled libs
wget https://github.com/SudoDEM/3rdlibs/releases/latest/download/3rdlibs-${OS}.tar.xz
tar xf 3rdlibs-${OS}.tar.xz

cd 3rdlib
#cp all libraries
cp -rf HeaderLib/lib/* 3rdlibs/
cp -r HeaderLib/lib/py/* 3rdlibs/py/
