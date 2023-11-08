#!/bin/bash

# system level dependencies
sudo apt-get install -y build-essential cmake zlib1g-dev freeglut3-dev libbz2-dev libglib2.0 python3-dev python3-pip python3-tk

# Qt5 based library
sudo apt-get install -y qt5-default libqt5svg5-dev libqt5webkit5-dev libqt5designer5

# python pip dependencies
python3 -m pip install numpy matplotlib xlib ipython


