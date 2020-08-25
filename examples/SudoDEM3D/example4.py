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
from sudodem._superquadrics_utils import *
import math
import random
import numpy as np

a=1.1734e-2
isSphere=False
ap=0.4
eps=0.5
num_s=200
num_t=8000
trials=0


wallboxsize=[10e-1,0,0] #size of the layer for the particles generation
boxsize=np.array([7e-1,7e-1,20e-1])

def down():
    for b in O.bodies:
        if(isinstance(b.shape,PolySuperellipsoid)):
            if(len(b.intrs())==0):
                b.state.vel=[0,0,-2]
                b.state.angVel=[0,0,0]

def gettop():
    zmax=0
    for b in O.bodies:
        if(isinstance(b.shape,PolySuperellipsoid)):
            if(b.state.pos[2]>=zmax):
                zmax=b.state.pos[2]
    return zmax


def Genparticles(a,eps,ap,mat,boxsize,num_s,num_t,bottom):
    global trials
    gap=max(a,a*ap)
    #print(gap)
    num=0
    coor=[]
    width=(boxsize[0]-2.*gap)/2.
    length=(boxsize[1]-2.*gap)/2.
    height=(boxsize[2]-2.*gap)
    #iteration=0

    while num<num_s and trials<num_t:
        isOK=True
        pos=[0]*3
        pos[0]=random.uniform(-width,width)
        pos[1]=random.uniform(gap-length,length-gap)
        pos[2]=random.uniform(bottom+gap,bottom+gap+height)

        for i in range(0,len(coor)):
            distance=sum([((coor[i][j]-pos[j])**2.) for j in range(0,3)])
            if(distance<(4.*gap*gap)):

                isOK=False
                break

        if(isOK==True):
            coor.append(pos)
            num+=1
            trials+=1
            #print(num)
    return coor


def Addlayer(a,eps,ap,mat,boxsize,num_s,num_t):
   bottom=gettop()
   coor=Genparticles(a,eps,ap,mat,boxsize,num_s,num_t,bottom)
   for b in coor:
      #xyz = np.array([1.0,0.5,0.8,1.5,1.2,0.4])*0.1
      xyz = np.array([1.0,0.5,0.8,1.5,0.2,0.4])*0.1
      bb=NewPolySuperellipsoid([1.0,1.0],xyz,mat,True,isSphere)#
      bb.state.pos=[b[0],b[1],b[2]]
      bb.state.vel=[0,0,-1]
      O.bodies.append(bb)
   down()



p_mat = PolySuperellipsoidMat(label="mat1",Kn=1e5,Ks=7e4,frictionAngle=math.atan(0.3),density=2650,betan=0,betas=0) #define Material with default values
wall_mat = PolySuperellipsoidMat(label="wallmat",Kn=1e6,Ks=7e5,frictionAngle=0.0,betan=0,betas=0) #define Material with default values
wallmat_b = PolySuperellipsoidMat(label="wallmat",Kn=1e6,Ks=7e5,frictionAngle=math.atan(1),betan=0,betas=0) #define Material with default values

O.materials.append(p_mat)
O.materials.append(wall_mat)
O.materials.append(wallmat_b)


O.bodies.append(utils.wall(-0.5*wallboxsize[0],axis=0,sense=1, material = wall_mat))#left x
O.bodies.append(utils.wall(0.5*wallboxsize[0],axis=0,sense=-1, material = wall_mat))#right x
O.bodies.append(utils.wall(-0.5*boxsize[1],axis=1,sense=1, material = wall_mat))#front y
O.bodies.append(utils.wall(0.5*boxsize[1],axis=1,sense=-1, material = wall_mat))#back y

O.bodies.append(utils.wall(0,axis=2,sense=1, material =wallmat_b))#bottom z
#O.bodies.append(utils.wall(boxsize[0],axis=2,sense=-1, material = wallmat_b))#top z

newton=NewtonIntegrator(damping = 0.3,gravity=(0.,0.,-9.8),label="newton",isSuperquadrics=2)#isSuperquadrics=2 for Poly-superellipsoids

O.engines=[
   ForceResetter(),
   InsertionSortCollider([Bo1_PolySuperellipsoid_Aabb(),Bo1_Wall_Aabb()],verletDist=0.2*0.1),
   InteractionLoop(
      [Ig2_Wall_PolySuperellipsoid_PolySuperellipsoidGeom(), Ig2_PolySuperellipsoid_PolySuperellipsoid_PolySuperellipsoidGeom()],
      [Ip2_PolySuperellipsoidMat_PolySuperellipsoidMat_PolySuperellipsoidPhys()], # collision "physics"
      [PolySuperellipsoidLaw()]   # contact law -- apply forces
   ),
   newton
]
O.dt=5e-5
#adding particles
Addlayer(a,eps,ap,p_mat,boxsize,num_s,num_t)#Note: these particles may intersect with each other when we generate them.

#setting the Display properties. You can go to the Display panel to set them directly.
sudodem.qt.Gl1_PolySuperellipsoid.wire=False
sudodem.qt.Gl1_PolySuperellipsoid.slices=10
sudodem.qt.Gl1_Wall.div=0 #hide the walls
