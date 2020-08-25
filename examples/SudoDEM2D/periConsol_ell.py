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
from sudodem._superellipse_utils import *
from minieigen import *
import numpy as np
import math, random

random.seed(11245642)
####################################################
#########constants
####################################################
pos_filename = 'ballinfo-large-poly.txt'
isSphere = False
size_alpha = 1.0#1.15325632818805
r = 0.005 #particle size
Particle_Num = 400
num_in_1d = int(math.sqrt(Particle_Num))
binsize = r*2
boxsize=[binsize*num_in_1d,binsize*num_in_1d*2]


#O.cell.hSize = Matrix2(boxsize[0],0, 0,boxsize[1])
O.periodic=True
#O.periodicCell=((-0.5*boxsize[0],0.5*boxsize[1]),(0.5*boxsize[0],0.5*boxsize[1]))  # activates periodic boundary
O.cell.setBox(boxsize)
#O.periodicCell=() # deactivates periodic boundary
#############################
mat=SuperellipseMat(label="mat1",Kn=1e7,Ks=1e7,frictionAngle=math.atan(0.5),density=2650)
O.materials.append(mat)

def setContactFriction(friction):
	O.materials[0].frictionAngle=math.atan(friction)
	for itr in O.interactions:
		itr.phys.tangensOfFrictionAngle = friction

def gen_particles(num,r_ball,binsize):
    for j in range(num*2):
        y = r_ball + j*binsize
        for i in range(num):
            x = r_ball + i*binsize
            rr = r_ball#*random.uniform(0.5,1.0)
            body = NewSuperellipse_rot(rr,rr*0.5,1.0,mat,0.0/180.0*np.pi, False,z_dim=0.075)#body = NewSuperellipse(r_ball,r_ball*0.5,1.0,mat,True,False)
            O.bodies.append(body)
            O.bodies[-1].state.pos=(x,y)

def print_fnt():
	for itr in O.bodies[399].intrs():
		fn = itr.phys.normalForce
		ft = itr.phys.shearForce
		mu = ft.norm()/fn.norm()
		print mu

gen_particles(num_in_1d,r,binsize)

def blockRotation():
    for b in O.bodies:
        b.state.blockedDOFs = 'Z'

def post_consol():
   O.engines = O.engines[:4]
   #print getStress()
   print O.cell.hSize
   setContactFriction(0.5)
   O.cell.trsf=Matrix2.Identity
   O.cell.velGrad=Matrix2.Zero
   for p in O.bodies:
      p.state.vel = Vector2.Zero
      p.state.angVel = 0.0
      p.state.refPos = p.state.pos
      p.state.refOri = p.state.ori
   O.save('con_ell.xml.bz2')
   O.pause()


#
#print Matrix2.Identity
PT= PeriTriaxController(
      dynCell=True,
      z_dim = 0.075,
      goal=(-0.95e5,-0.95e5),
      stressMask=3,
      relStressTol=.001,
      maxUnbalanced=.001,
      maxStrainRate=(.5,.5),
      doneHook='post_consol1()',
      label='biax'
   )

PT2= PeriTriaxController(
      dynCell=True,
      z_dim = 0.075,
      goal=(-1e5,-1e5),
      stressMask=3,
      relStressTol=.001,
      maxUnbalanced=.001,
      maxStrainRate=(.5,.5),
      doneHook='post_consol()',
      label='biax'
   )

def post_consol1():
    print "State one over!"
    O.engines[4].doneHook='post_consol()'
    O.engines[4].goal=(-1e5,-1e5)
    for b in O.bodies:
        b.state.blockedDOFs=''
    print "FTPO=",utils.getFabricTensorPO2D()
    print "FTCN=",utils.getFabricTensorCN2D()
    print "MCN=",utils.getMeanCN()
    #O.pause()


blockRotation()

O.dt = 1e-5
newton=NewtonIntegrator(damping = 0.2,gravity=(0.,0.0),label="newton",isSuperellipse=True)
O.engines=[
   ForceResetter(),
   InsertionSortCollider([Bo1_Superellipse_Aabb()],verletDist=1e-3),
   InteractionLoop(
	  [Ig2_Superellipse_Superellipse_SuperellipseGeom()],
	  [Ip2_SuperellipseMat_SuperellipseMat_SuperellipsePhys()], # collision "physics"
	  [SuperellipseLaw()]   # contact law -- apply forces
   ),
   newton,
   PT
]
