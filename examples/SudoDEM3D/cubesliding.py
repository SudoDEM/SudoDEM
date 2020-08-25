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
#This is a script for testing sliding of a cube on a wall
#Note that we can not change the decline angle of the wall rather than a facet
#But we prefer to use a wall here
#Hence, as a workaround, we just decompose the gravitational acceleration into x and z axises. 
from sudodem import  _gjkparticle_utils
from sudodem import *

from sudodem._gjkparticle_utils import *
import random as rand
import math

cube_v0 = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
cube_v = [[i*0.01 for i in j] for j in cube_v0]

########define some parameters#################
isSphere=False

inclineAngle = 30.0/180.0*math.pi
R = 0.01

fname = 'state.dat'

def outputData(fname, clearFile = False):
    '''
    output the position and velocity of the particle
    '''
    b = O.bodies[-1].state
    f = None
    if(clearFile):
        f = open(fname, 'w')
        print>>f, "iteration, position, velocity"
    else:
        f = open(fname, 'a')
    print>>f, O.iter, b.pos[0], b.vel[0]
    f.close()

p_mat = GJKParticleMat(label="mat1",Kn=1e5,Ks=1e5,frictionAngle=math.atan(1.0),density=2650,betan=0,betas=0) #define Material with default values
wallmat_b = GJKParticleMat(label="wallmat",Kn=1e5,Ks=1e5,frictionAngle=math.atan(1.0),betan=0.5,betas=0.5) #define Material with default values

O.materials.append(p_mat)
O.materials.append(wallmat_b)

# create the box and add it to O.bodies
O.bodies.append(utils.wall(0,axis=2,sense=1, material =wallmat_b))#bottom wall along z axis

body = GJKCuboid([0.005,0.005,0.005],0.05*1e-2,p_mat,False)
gra = 1e-6*2650*9.8*math.cos(inclineAngle)/1e5
body.state.pos=[0.0,0.0,0.005-gra]
O.bodies.append(body)
O.bodies[-1].shape.color=(rand.random(),rand.random(),rand.random())


newton=NewtonIntegrator(damping = 0.0,gravity=(9.8*math.sin(inclineAngle),0.0,-9.8*math.cos(inclineAngle)),label="newton",isSuperquadrics=4)


O.engines=[
   ForceResetter(),
   InsertionSortCollider([Bo1_GJKParticle_Aabb(),Bo1_Wall_Aabb()],verletDist=0.2*0.01),
   InteractionLoop(
      [Ig2_Wall_GJKParticle_GJKParticleGeom(), Ig2_GJKParticle_GJKParticle_GJKParticleGeom()],
      [Ip2_GJKParticleMat_GJKParticleMat_GJKParticlePhys()], # collision "physics"
      [GJKParticleLaw()]   # contact law -- apply forces
   ),
   newton,
   PyRunner(command='outputData(fname)',iterPeriod=100,dead = False)
]
O.dt=1e-6

outputData(fname, clearFile = True)
