# HOW TO USE THE SCRIPTS FOR QUICK INSTALLATION and COMPILATION
## Using compiled executables
You can download the compiled executables from the [release page](https://github.com/SudoDEM/3rdlibs/releases/latest) on GitHub. Or type
the following commands in a terminal to download them.
```bash
wget https://github.com/SudoDEM/3rdlibs/releases/latest/download/SudoDEM2D-ubuntu-20.04.tar.xz
wget https://github.com/SudoDEM/3rdlibs/releases/latest/download/SudoDEM3D-ubuntu-20.04.tar.xz
```
## Using compiled 3rdlibs
We have compiled the 3rd-party libraries on Ubuntu 20.04 so that you can reuse them.
1. Install some tools for compilation, e.g., gcc, cmake, etc. You will be asked for the root permission.
```bash
./install3rd_base.sh
```
2. Download the compiled 3rd-party libraries
```bash
./install3rdcore_ubuntu20.04.sh
```
3. Compile SudoDEM2D
```bash
./firstCompile_SudoDEM2D.sh
```
4. Compile SudoDEM3D
```bash
./fristCompile_SudoDEM3D.sh
```
## Compiling everything from sources
You may need to compile the 3rd-party libraries for other Linux distributions/versions.
1. Install some tools for compilation, e.g., gcc, cmake, etc. You will be asked for the root permission.
```bash
./install3rd_base.sh
```
2. Compile the 3rd-party libraries
```bash
./install3rdcore_sudosimlab.sh
```
3. Compile SudoDEM2D and SudoDEM3D
```bash
./firstCompile_SudoDEM2D.sh
./fristCompile_SudoDEM3D.sh
```
## Compiling everything from sources (sudosimlab)
You may encounter slow downloading of the sources from GitHub. Please swtich to our SudoSimLab server. Refer to "install3rdcore_sudosimlabServer.sh".

