#!/bin/bash
cd ..
WORKSPACE=$PWD
#6. Compile SudoDEM3D
mkdir build3d
cd build3d
cmake ../SudoDEM3D
#make
make -j3 #use 3 threads to compile
make install
#the compiled files will be install at `sudodeminstall/SudoDEM3D/'.
#copy the dynamic files `3rdlibs' into `sudodeminstall/SudoDEM3D/lib/'
cd $WORKSPACE
cp -rf 3rdlib/3rdlibs sudodeminstall/SudoDEM3D/lib/
#make sure to copy IPython pakage into `sudodeminstall/SudoDEM3D/lib/py'.

#7. Set PATH in ~/.bashrc
#Append the following command into `~/.bashrc'.
#PATH=${PATH}:/home/xxx/SudoDEM/sudodeminstall/SudoDEM3D/bin
echo 'PATH=${PATH}:'$WORKSPACE/sudodeminstall/SudoDEM3D/bin >> ~/.bashrc
