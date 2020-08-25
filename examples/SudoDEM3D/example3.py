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
from sudodem import plot, _gjkparticle_utils
from sudodem import *

from sudodem._gjkparticle_utils import *
import random as rand
import math

cube_v0 = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
cube_v = [[i*0.01 for i in j] for j in cube_v0]

########define some parameters#################
isSphere=False

num_x = 5
num_y = 5
num_z = 20
R = 0.01

#######define some auxiliary functions#########

#define a latice grid
#R: distance between two neighboring nodes
#num_x: number of nodes along x axis
#num_y: number of nodes along y axis
#num_z: number of nodes along z axis
#return a list of the postions of all nodes
def GridInitial(R,num_x=10,num_y=10,num_z=20):
   pos = list()
   for i in range(num_x):
      for j in range(num_y):
         for k in range(num_z):
            x = i*R*2.0
            y = j*R*2.0
            z = k*R*2.0
            pos.append([x,y,z])
   return pos

#generate a sample
def GenSample(r,pos):
   for p in pos:

		#body = GJKPolyhedron([],[1.e-2,1.e-2,1.e-2],0.05*1e-2,p_mat,False)
		body = GJKCuboid([0.005,0.005,0.005],0.05*1e-2,p_mat,True)
		#body = GJKCone(0.005,0.01,0.05*1e-2,p_mat,True)
		#body = GJKCylinder(0.005,0.01,0.05*1e-2,p_mat,True)
		body.state.pos=p
		O.bodies.append(body)
		O.bodies[-1].shape.color=(rand.random(),rand.random(),rand.random())


p_mat = GJKParticleMat(label="mat1",Kn=1e4,Ks=7e3,frictionAngle=math.atan(0.5),density=2650,betan=0,betas=0) #define Material with default values
#mat = GJKParticleMat(label="mat1",Kn=1e6,Ks=1e5,frictionAngle=0.5,density=2650) #define Material with default values
wall_mat = GJKParticleMat(label="wallmat",Kn=1e4,Ks=7e3,frictionAngle=0.0,betan=0,betas=0) #define Material with default values
wallmat_b = GJKParticleMat(label="wallmat",Kn=1e4,Ks=7e3,frictionAngle=math.atan(1),betan=0,betas=0) #define Material with default values

O.materials.append(p_mat)
O.materials.append(wall_mat)
O.materials.append(wallmat_b)

# create the box and add it to O.bodies
O.bodies.append(utils.wall(-R,axis=0,sense=1, material = wall_mat))#left wall along x axis
O.bodies.append(utils.wall(2.0*R*num_x-R,axis=0,sense=-1, material = wall_mat))#right wall along x axis
O.bodies.append(utils.wall(-R,axis=1,sense=1, material = wall_mat))#front wall along y axis
O.bodies.append(utils.wall(2.0*R*num_y-R,axis=1,sense=-1, material = wall_mat))#back wall along y axis
O.bodies.append(utils.wall(-R,axis=2,sense=1, material =wallmat_b))#bottom wall along z axis


# create a lattice grid
pos = GridInitial(R,num_x=num_x,num_y=num_y,num_z=num_z) # get positions of all nodes
# create particles at each nodes
GenSample(R,pos)


newton=NewtonIntegrator(damping = 0.3,gravity=(0.,0.,-9.8),label="newton",isSuperquadrics=4)

O.engines=[
   ForceResetter(),
   InsertionSortCollider([Bo1_GJKParticle_Aabb(),Bo1_Wall_Aabb()],verletDist=0.2*0.01),
   InteractionLoop(
      [Ig2_Wall_GJKParticle_GJKParticleGeom(), Ig2_GJKParticle_GJKParticle_GJKParticleGeom()],
      [Ip2_GJKParticleMat_GJKParticleMat_GJKParticlePhys()], # collision "physics"
      [GJKParticleLaw()]   # contact law -- apply forces
   ),
   newton
]
O.dt=1e-5
