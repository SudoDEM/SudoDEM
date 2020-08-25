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
from sudodem import utils
from minieigen import *
import numpy as np
import math

loading_stress = [100]#
saving_num = 0
#############################
####################################################
def setContactFriction(friction):
	O.materials[0].frictionAngle=math.atan(friction)
	for itr in O.interactions:
		itr.phys.tangensOfFrictionAngle = friction
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
O.load("con_disk"+suffix+".xml.bz2")


O.cell.velGrad=Matrix2(0,1,0,0)#distortion along y, simple shear
newton=NewtonIntegrator(damping = 0.2,gravity=(0.,0.),label="newton")
O.engines=[
   ForceResetter(),
   InsertionSortCollider([Bo1_Disk_Aabb(),Bo1_Wall_Aabb()],verletDist=0.25e-3),
   InteractionLoop(
	  [Ig2_Wall_Disk_ScGeom(), Ig2_Disk_Disk_ScGeom()],
	  [Ip2_FrictMat_FrictMat_FrictPhys()], # collision "physics"
	  [Law2_ScGeom_FrictPhys_CundallStrack()]   # contact law -- apply forces
   ),
   newton,
   #PT,
   #PyRunner(command='output()',iterPeriod=3000,label='outdata',dead = False),
   #PyRunner(command='savedata()',realPeriod=3600,label='savedata',dead = False)#saving data every two hours
]



O.dt = 1e-4

setContactFriction(0.5)
#O.cell.trsf=Matrix2.Identity
#O.cell.velGrad=Matrix2.Zero
