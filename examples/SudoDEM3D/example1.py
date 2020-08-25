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
###############################################
# granular packing of super-ellipsoids
# We position particles at a lattice grid,
# then particles free fall into a cubic box.
###############################################
# import some modules
from sudodem import _superquadrics_utils

from sudodem._superquadrics_utils import *
import math
import random as rand
import numpy as np

########define some parameters#################
isSphere=False

num_x = 5
num_y = 5
num_z = 20
R = 0.1

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
		epsilon1 = rand.uniform(0.7,1.6)
		epsilon2 = rand.uniform(0.7,1.6)
		a = r
		b = r*rand.uniform(0.4,0.9)#aspect ratio
		c = r*rand.uniform(0.4,0.9)

		body = NewSuperquadrics2(a,b,c,epsilon1,epsilon2,p_mat,True,isSphere)
		body.state.pos=p
		O.bodies.append(body)

#########setup a simulation####################
# material for particles
p_mat = SuperquadricsMat(label="mat1",Kn=1e5,Ks=7e4,frictionAngle=math.atan(0.3),density=2650,betan=0,betas=0)
# betan and betas are coefficients of viscous damping at contact, no viscous damping with 0 by default.
# material for the side walls
wall_mat = SuperquadricsMat(label="wallmat",Kn=1e6,Ks=7e5,frictionAngle=0.0,betan=0,betas=0)
# material for the bottom wall
wallmat_b = SuperquadricsMat(label="wallmat",Kn=1e6,Ks=7e5,frictionAngle=math.atan(1),betan=0,betas=0)
# add materials to O.materials
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

# create engines
# define a Newton engine
newton=NewtonIntegrator(damping = 0.1,gravity=(0.,0.,-9.8),label="newton",isSuperquadrics=1) # isSuperquadrics: 1 for superquadrics

O.engines=[
   ForceResetter(), InsertionSortCollider([Bo1_Superquadrics_Aabb(),Bo1_Wall_Aabb()],verletDist=0.2*0.1),
   InteractionLoop(
      [Ig2_Wall_Superquadrics_SuperquadricsGeom(),    Ig2_Superquadrics_Superquadrics_SuperquadricsGeom()],
      [Ip2_SuperquadricsMat_SuperquadricsMat_SuperquadricsPhys()], # collision "physics"
      [SuperquadricsLaw()]   # contact law
   ),
   newton

]
O.dt=5e-5
