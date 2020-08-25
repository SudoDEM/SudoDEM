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
import numpy as np
import math
import random
random.seed(11245642)
####################################################
#########constants
####################################################
isSphere = False
size_alpha = 1.0#1.15325632818805
boxsize=[50.0,10.0]
r = 0.5 #particle size
num_s=100
num_t=8000
trials=0
boxsize=[50.0,10.0]

#boxsize=[2.4,2.4,2.4]
#loading_stress = [50,100,250]#
loading_stress = [100]#
#############################
mat =SuperellipseMat(label="mat1",Kn=1e8,Ks=7e7,frictionAngle=math.atan(0.225),density=2650)
O.materials.append(mat)

O.bodies.append(utils.wall(-0.5*boxsize[1],axis=1,sense=1, material = 'mat1'))#ground
O.bodies.append(utils.wall(-0.5*boxsize[0],axis=0,sense=1, material = 'mat1'))#left
O.bodies.append(utils.wall(0.5*boxsize[0],axis=0,sense=-1, material = 'mat1'))#right

###
def gettop():
    ymax=0
    for b in O.bodies:
        if(isinstance(b.shape,Superellipse)):
            if(b.state.pos[1]>=ymax):
                ymax=b.state.pos[1]
    return ymax


def Genparticles(r,mat,boxsize,num_s,num_t,bottom):
    global trials
    gap=r
    #print(gap)
    num=0
    coor=[]
    width=(boxsize[0]-2.*gap)/2.
    height=(boxsize[1]-2.*gap)
    #iteration=0

    while num<num_s and trials<num_t:
        isOK=True
        pos=[0]*2
        pos[0]=random.uniform(-width,width)
        pos[1]=random.uniform(bottom+gap,bottom+gap+height)

        for i in range(0,len(coor)):
            distance=sum([((coor[i][j]-pos[j])**2.) for j in range(0,2)])
            if(distance<(4.*gap*gap)):

                isOK=False
                break

        if(isOK==True):
            coor.append(pos)
            num+=1
            trials+=1
            #print(num)
    return coor


def Addlayer(r,mat,boxsize,num_s,num_t):
    bottom=gettop()
    coor=Genparticles(r,mat,boxsize,num_s,num_t,bottom)
    for b in coor:
        #bb = utils.sphere(b,0.5,material = 'mat1')
	rr = r*random.uniform(0.5,1.0)
	bb = NewSuperellipse(1.5*r,r,1.0,mat,True,False)
	bb.state.pos=b
        bb.state.vel=[0,-1]
        O.bodies.append(bb)
	if len(O.bodies) - 3 > 2000:
		O.engines[-1].dead=True

O.dt = 1e-3

def out_energy():
	v1 = O.bodies[-1].state.vel
	v2 = O.bodies[-2].state.vel
	print v1.norm() + v2.norm()
newton=NewtonIntegrator(damping = 0.1,gravity=(0.,-10.0),label="newton",isSuperellipse=True)

O.engines=[
   ForceResetter(),
   InsertionSortCollider([Bo1_Superellipse_Aabb(),Bo1_Wall_Aabb()],verletDist=0.1),
   InteractionLoop(
      [Ig2_Wall_Superellipse_SuperellipseGeom(), Ig2_Superellipse_Superellipse_SuperellipseGeom()],
      [Ip2_SuperellipseMat_SuperellipseMat_SuperellipsePhys()], # collision "physics"
      [SuperellipseLaw()]   # contact law -- apply forces
   ),
   #GravityEngine(gravity=(0,0,-9.81)),
   #PyRunner(command='setvel1()',iterPeriod=1000,label='setveler1',dead=False),
   #triax,
   newton,
   PyRunner(command='Addlayer(r,mat,boxsize,num_s,num_t)',virtPeriod=0.1,label='check',dead = False)
]
