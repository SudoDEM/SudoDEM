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
import random as rand
rand.seed(11245642)
####################################################
epsilon = 1.4
isSphere = False
boxsize=[16.2e-3,16.2e-3,16.2e-3]
############################
#some auxiliary functions
#generate a sample of monodisperse superellipsoids with random orientaions and positions
def gen_sample(width, depth, height, r_max = 1.0e-3):
   b_num_tot = 5000
   for j in range(b_num_tot):
      rb = r_max*0.5
      x = rand.uniform(rb, width-rb)
      y = rand.uniform(rb, depth - rb)
      z = rand.uniform(rb, height - rb)
      r = r_max*0.5
      b = NewSuperquadrics2(r*1.25,r,r,epsilon,epsilon,mat,True,isSphere)
      b.state.pos=(x, y, z)
      O.bodies.append(b)
#Materials definition
mat = SuperquadricsMat(label="mat1",Kn=3e4,Ks=3e4,frictionAngle=math.atan(0.1),density=2650e4) #define Material with default values
wallmat1 = SuperquadricsMat(label="wallmat1",Kn=1e6,Ks=0.,frictionAngle=0.) #define Material with default values
wallmat2 = SuperquadricsMat(label="wallmat2",Kn=1e6,Ks=0.,frictionAngle=0.) #define Material with default values
O.materials.append(mat)
O.materials.append(wallmat1)
O.materials.append(wallmat2)
#############################################
#####create confining walls
############################################
O.bodies.append(utils.wall(0,axis=0,sense=1, material = 'wallmat1'))#left x
O.bodies.append(utils.wall(boxsize[0],axis=0,sense=-1, material = 'wallmat1'))#right x
O.bodies.append(utils.wall(0,axis=1,sense=1, material = 'wallmat1'))#front y
O.bodies.append(utils.wall(boxsize[1],axis=1,sense=-1, material = 'wallmat1'))#back y

O.bodies.append(utils.wall(0,axis=2,sense=1, material = 'wallmat2'))#bottom z
O.bodies.append(utils.wall(boxsize[2],axis=2,sense=-1, material = 'wallmat2'))#top z
#####create particles

gen_sample(boxsize[0], boxsize[1], boxsize[2])

####function for saving data (e.g., simulation states here)
def savedata():
        O.save(str(O.time)+'.xml.bz2')

#########################################
####triax engine
#CAUTION: the CompressionEngine regards the first six bodies as confining walls by default. So the walls should be generated first.
triax=CompressionEngine(
   goalx = 1.0e5, # confining stress 100 kPa
   goaly = 1.0e5,
   goalz = 1.0e5,
	savedata_interval = 10000,
	echo_interval = 2000,
	continueFlag = False,
	max_vel = 10.0,
   gain_alpha = 0.5,
   f_threshold = 0.01, #1-f/goal
   fmin_threshold = 0.01, #the ratio of deviatoric stress to mean stress
   unbf_tol = 0.01, #unblanced force ratio
	file='record-con',
)

############################



newton=NewtonIntegrator(damping = 0.3,gravity=(0.,0.,0.),label="newton",isSuperquadrics=1)
O.engines=[
   ForceResetter(),
   InsertionSortCollider([Bo1_Superquadrics_Aabb(),Bo1_Wall_Aabb()],verletDist=0.1e-3),
   InteractionLoop(
      [Ig2_Wall_Superquadrics_SuperquadricsGeom(), Ig2_Superquadrics_Superquadrics_SuperquadricsGeom()],
      [Ip2_SuperquadricsMat_SuperquadricsMat_SuperquadricsPhys()],
      [SuperquadricsLaw()]
   ),
   triax,
   newton,
   PyRunner(command='quiet_ball()',iterPeriod=2,iterLast = 50000,label='calm',dead = False),
   PyRunner(command='savedata()',realPeriod=7200.0,label='savedata',dead = False)#saving data every two hours
]


#delete the balls outside the box
def del_balls():
        num_del = 0
        #wall position
        left = O.bodies[0].state.pos[0]#x
        right = O.bodies[1].state.pos[0]
        front = O.bodies[2].state.pos[1]#y
        back = O.bodies[3].state.pos[1]
        bottom = O.bodies[4].state.pos[2]#z
        top = O.bodies[5].state.pos[2]
        for b in O.bodies:
                if isinstance(b.shape,Superquadrics):
                        (x,y,z) = b.state.pos
                        if x > right or x < left or y < front or y > back or z < bottom or z > top:
                                O.bodies.erase(b.id)
                                num_del += 1
        print str(num_del)+" particles are deleted!"



def quiet_ball():
    global calm_num

    newton.quiet_system_flag = True
    if calm_num > 2000:
        O.engines[-2].iterPeriod= 5
        calm_num += 5
    else:
        calm_num += 2
    if calm_num > 40000:
        if calm_num < 40010:
            print "calm procedure is over"
            del_balls()
        if calm_num > 50000:
            O.engines[-2].dead = True
            O.save("init_assembly.xml.bz2")
            #consolidation begins
            triax.dead=False
            print "consolidation begins"


O.dt=1e-4


triax.dead=True
calm_num = 0
