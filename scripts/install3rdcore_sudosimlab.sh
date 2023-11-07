#!/bin/bash

# python pip dependencies
python3 -m pip install numpy matplotlib xlib ipython

# 1 create the 3rdlib dir

cd ..
mkdir 3rdlib

cd 3rdlib

mkdir 3rdlibs
mkdir -p 3rdlibs/py


# 2 workspace = sudodem_3rdlib
WORKSPACE=$PWD

# 3 Boost install
#wget https://www.sudosimlab.com/downloadData/3rdlibs/boost_1_67_0.tar.gz
#tar xzf boost_1_67_0.tar.gz
wget https://boostorg.jfrog.io/artifactory/main/release/1.67.0/source/boost_1_67_0.tar.gz
tar xzf boost_1_67_0.tar.gz
cd boost_1_67_0
./bootstrap.sh  --with-libraries=python,thread,filesystem,iostreams,regex,serialization,system,date_time link=shared runtime-link=shared --without-icu --with-python=python3
./b2 −j4 #compile with 3 threads
./b2 install --prefix=$WORKSPACE/HeaderLib

# 4 Eigen and minieigen
cd $WORKSPACE
wget https://www.sudosimlab.com/downloadData/3rdlibs/eigen-3.3.5.tar.gz
tar xzf eigen-3.3.5.tar.gz
mv eigen-git-mirror-3.3.5 Eigen-3.3.5

wget https://www.sudosimlab.com/downloadData/3rdlibs/minieigen_1_py3.tar.xz
tar xf minieigen_1_py3.tar.xz
cd minieigen_1_py3
mkdir build
cd build 
cmake .. -DCMAKE_INSTALL_PREFIX=$WORKSPACE/HeaderLib/lib/py
make -j4
make install

# 5 LibQGLViewer-2.6.3
cd $WORKSPACE
wget https://www.sudosimlab.com/downloadData/3rdlibs/libQGLViewer-2.6.3.tar.gz
tar xzf libQGLViewer-2.6.3.tar.gz
cd libQGLViewer-2.6.3-2.6.3/QGLViewer
qmake 
make
make install


# 6 Sip
cd $WORKSPACE
wget https://www.sudosimlab.com/downloadData/3rdlibs/sip_4_19_13.tar.xz
tar xf sip_4_19_13.tar.xz
cd sip_4_19_13
mkdir -p $WORKSPACE/../SudoDEMTools/sip/sip2pyqt
python3 configure.py --sip-module=PyQt5.sip --bindir=$WORKSPACE/../SudoDEMTools/sip/bin --destdir=$WORKSPACE/HeaderLib/lib/py --incdir=$WORKSPACE/../SudoDEMTools/sip/include --sipdir=$WORKSPACE/../SudoDEMTools/sip/sip2pyqt
make -j4
make install

# 7 PyQt5
cd $WORKSPACE
wget https://www.sudosimlab.com/downloadData/3rdlibs/pyqt_5_11_3.tar.xz
tar xf pyqt_5_11_3.tar.xz
cd pyqt_5_11_3

python3 configure.py --qmake=/usr/bin/qmake --sip=$WORKSPACE/../SudoDEMTools/sip/bin/sip --sip-incdir=$WORKSPACE/../SudoDEMTools/sip/include --bindir=$WORKSPACE/../SudoDEMTools/PyQt5/bin --designer-plugindir=$WORKSPACE/../SudoDEMTools/PyQt5/plugin --sipdir=$WORKSPACE/../SudoDEMTools/sip/sip2pyqt --destdir=$WORKSPACE/HeaderLib/lib/py --stubsdir=$WORKSPACE/../SudoDEMTools/PyQt5/stubs --qml-plugindir=$WORKSPACE/../SudoDEMTools/PyQt5/qml --confirm-license
make -j4
make install


# 8 cp all libraries
cd $WORKSPACE
cp -rf HeaderLib/lib/* 3rdlibs/
cp -r HeaderLib/lib/py/* 3rdlibs/py/
