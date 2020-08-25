"""
 #     SudoDEM : A Discrete Element Code for Non-spherical Particles
 #
 #               Copyright (C) 2016 - 2020   Shiwei Zhao
 #
 #
 #		      South China University of Technology
 #       Hong Kong University of Science and Technology
 #
 #  This script is free software; you can redistribute it and/or
 #  modify it under the terms of the GNU General Public
 #  License as published by the Free Software Foundation; either
 #  version 3.0 of the License, or (at your option) any later version.
 #
 #  This program is distributed in the hope that it will be useful,
 #  but WITHOUT ANY WARRANTY; without even the implied warranty of
 #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 #  General Public License for more details.
 #
 #  You should have received a copy of the GNU General Public
 #  License along with this library; if not, write to the Free Software
 #  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 #
 """
##################################################
from sudodem import plot, _superquadrics_utils

from sudodem._superquadrics_utils import *

import math

saving_num = 0
####################################################
def savedata():
        global saving_num
        saving_num += 1
        #print "testing, no data saved"
        O.save('shear'+str(saving_num)+'.xml.bz2')

def quiet_system():
    for i in O.bodies:
        i.state.vel=(0.,0.,0.)
        i.state.angVel=(0.,0.,0.)
suffix=""
O.load("final_con"+suffix+".xml.bz2")
triax = O.engines[3]
O.dt = 5e-5
friction = 0.5
O.materials[0].frictionAngle=math.atan(friction)
if 1:
   for itr in O.interactions:
      id = min(itr.id1,itr.id2)
      if id > 5:
         itr.phys.tangensOfFrictionAngle = friction
         #itr.phys.betan = 0.5

#########################
#shear begins
#########################
def shear():
    quiet_system()
    triax.z_servo = False
    triax.ramp_interval = 10
    triax.ramp_chunks = 2
    triax.file="sheardata"+suffix
    triax.goalz = 0.01
    triax.goalx = 1e5
    triax.goaly = 1e5
    triax.savedata_interval = 2500
    triax.echo_interval = 2000
    triax.target_strain = 0.4

shear()

O.engines[-1].iterPeriod = 2500
O.engines[-1].realPeriod = 0.0
