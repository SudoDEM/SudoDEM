#!/bin/bash
cd ..
WORKSPACE=$PWD
#6. Compile SudoDEM2D
mkdir build2d
cd build2d
cmake ../SudoDEM2D
#make
make -j3 #use 3 threads to compile
make install
#the compiled files will be install at `sudodeminstall/SudoDEM2D/'.
#copy the dynamic files `3rdlibs' into `sudodeminstall/SudoDEM2D/lib/'
cd $WORKSPACE
cp -rf 3rdlib/3rdlibs sudodeminstall/SudoDEM2D/lib/
#make sure to copy IPython pakage into `sudodeminstall/SudoDEM2D/lib/py'.

#7. Set PATH in ~/.bashrc
#Append the following command into `~/.bashrc'.
#PATH=${PATH}:/home/xxx/SudoDEM/sudodeminstall/SudoDEM2D/bin
echo 'PATH=${PATH}:'$WORKSPACE/sudodeminstall/SudoDEM2D/bin >> ~/.bashrc
