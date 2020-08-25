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
import math, random

random.seed(11245642)
####################################################
#########constants
####################################################
pos_filename = 'ballinfo-large-poly.txt'
epsilon = 1.0
isDisk = True
size_alpha = 1.0#1.15325632818805
boxsize=[0.5,0.5]
r = 0.5 #particle size
num_s=100
num_t=500
trials=0

Particle_Num = 400


GSD1 = {
0.5:	0.007811523,
0.595:	0.0175761999,
0.649:	0.0573008582,
0.672:	0.0670342374,
0.685:	0.0632896622,
0.695:	0.0605968541,
0.71:	0.1044659969,
0.725:	0.115410723,
0.737:	0.1016323632,
0.75:	0.0964385557,
0.765:	0.0908762049,
0.784:	0.0844279608,
0.796:	0.0403334423,
0.813:	0.0378558372,
0.85:	0.0369668044,
1.0:	0.0179827769
}

GSD = {
1.0:1.0
}

#O.cell.hSize = Matrix2(boxsize[0],0, 0,boxsize[1])
O.periodic=True
#O.periodicCell=((-0.5*boxsize[0],0.5*boxsize[1]),(0.5*boxsize[0],0.5*boxsize[1]))  # activates periodic boundary
O.cell.setBox(boxsize)
#O.periodicCell=() # deactivates periodic boundary
#############################
mat=FrictMat(label="mat1",Kn=1e6,Ks=1e6,frictionAngle=math.atan(0.001),density=2650)
O.materials.append(mat)
matwall=FrictMat(label="matwall",Kn=1e6,Ks=1e6,frictionAngle=math.atan(0.3),density=2650)
O.materials.append(matwall)

########walls
#O.bodies.append(utils.wall(-0.5*boxsize[1],axis=1,sense=1, material = 'matwall'))#ground
#O.bodies.append(utils.wall(-0.5*boxsize[0],axis=0,sense=1, material = 'matwall'))#left
#O.bodies.append(utils.wall(0.5*boxsize[0],axis=0,sense=-1, material = 'matwall'))#right
#O.bodies.append(utils.wall(0.5*boxsize[1],axis=1,sense=-1, material = 'matwall'))#ground
############
def setContactFriction(friction):
	O.materials[0].frictionAngle=math.atan(friction)
	for itr in O.interactions:
		itr.phys.tangensOfFrictionAngle = friction

def gettop():
    ymax=0
    for b in O.bodies:
        if(isinstance(b.shape,Disk)):
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
        bb = utils.disk(b,0.5,material = 'mat1')
        bb.state.vel=[0,-1]
        O.bodies.append(bb)

def gen_particles(b_num,r_ball, boxsize):
	for j in range(b_num):
		x = random.uniform(r_ball-0.5*boxsize[0], 0.5*boxsize[0]-r_ball)
		y = random.uniform(r_ball-0.5*boxsize[1], 0.5*boxsize[1] - r_ball)
		body = utils.disk((x,y),r_ball,material = 'mat1')
		O.bodies.append(body)

def gen_sample(boxsize, r_max = 1.0):
	b_num_tot = Particle_Num
	r = r_max
	r0 = r*0.5
	for i in GSD.keys():
		r_ball = i*0.5e-2
		b_num = int(b_num_tot*GSD[i])
		gen_particles(b_num, r_ball, boxsize)

def print_fnt():
	for itr in O.bodies[399].intrs():
		fn = itr.phys.normalForce
		ft = itr.phys.shearForce
		mu = ft.norm()/fn.norm()
		print mu

#particles
gen_sample(boxsize)

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
   O.save('con_disk.xml.bz2')
   O.pause()
#
#print Matrix2.Identity
PT= PeriTriaxController(
      dynCell=True,
      z_dim = 0.075,#dimension along z
      goal=(-1.e5,-1.e5),
      stressMask=3,
      relStressTol=.001,
      maxUnbalanced=.001,
      maxStrainRate=(.5,.5),
      doneHook='post_consol()',
      label='biax'
   )

O.dt = 1e-4
newton=NewtonIntegrator(damping = 0.2,gravity=(0.,0.0),label="newton")
O.engines=[
   ForceResetter(),
   InsertionSortCollider([Bo1_Disk_Aabb(),Bo1_Wall_Aabb()],verletDist=0.25e-2),
   InteractionLoop(
	  [Ig2_Wall_Disk_ScGeom(), Ig2_Disk_Disk_ScGeom()],
	  [Ip2_FrictMat_FrictMat_FrictPhys()], # collision "physics"
	  [Law2_ScGeom_FrictPhys_CundallStrack()]   # contact law -- apply forces
   ),
   newton,
   PT,
#PyRunner(command='Addlayer(r,mat,boxsize,num_s,num_t)',virtPeriod=0.1,label='check',dead = False)
]
