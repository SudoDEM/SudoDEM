# encoding: utf-8
# 2010 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>

'''
Motion of a "sinusoidal" beam made of cylinders. The test checks the position and velocity of the free end of the bending beam subjected to gravitational load. It is similar to scripts/test/chained-cylinder-spring.py but with less elements. positions and velocity are compared during the transient oscillations, only positions are compared for the larger time since residual velocity is compiler-dependent (see https://lists.launchpad.net/sudodem-dev/msg06178.html).
'''

import unittest
import random
from sudodem.wrapper import *
from sudodem._customConverters import *
from sudodem import utils
from math import *

try:
	from minieigen import *
except ImportError:
	from miniEigen import *

class TestCohesiveChain(unittest.TestCase):
	# prefix test names with PBC:
	def setUp(self):
		O.reset();
		young=1.0e3
		poisson=5
		density=2.60e3
		frictionAngle=radians(30)
		O.materials.append(CohFrictMat(young=young,poisson=poisson,density=density,frictionAngle=frictionAngle,normalCohesion=1e13,shearCohesion=1e13,momentRotationLaw=True))
		O.dt=1e-3
		O.engines=[
			ForceResetter(),
			InsertionSortCollider([
			Bo1_ChainedCylinder_Aabb(),
			Bo1_Disk_Aabb()]),
		InteractionLoop(
			[Ig2_ChainedCylinder_ChainedCylinder_ScGeom6D(),Ig2_Disk_ChainedCylinder_CylScGeom()],
			[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=True)],
			[Law2_ScGeom6D_CohFrictPhys_CohesionMoment()]),
		## Apply gravity
		## Motion equation
		NewtonIntegrator(damping=0.15,gravity=[0,-9.81,0])
		]
		#Generate a spiral
		Ne=10
		for i in range(0, Ne):
			omeg=95.0/float(Ne); hy=0.05; hz=0.07;
			px=float(i)*(omeg/60.0); py=sin(float(i)*omeg)*hy; pz=cos(float(i)*omeg)*hz;
			px2=float(i+1.)*(omeg/60.0); py2=sin(float(i+1.)*omeg)*hy; pz2=cos(float(i+1.)*omeg)*hz;
			utils.chainedCylinder(begin=Vector3(pz,py,px), radius=0.005,end=Vector3(pz2,py2,px2),color=Vector3(0.6,0.5,0.5))
		O.bodies[Ne-1].state.blockedDOFs='xyzXYZ'
	def testMotion(self):
		"CohesiveChain: velocity/positions tested in transient dynamics and equilibrium state"
		#target values
		tv1=-0.790576599652;tp1=-0.0483370740415;tv2=-0.000494271551993;tp2=-0.0611987315415;
		tolerance = 10e-3
		O.run(100,True)
		v1=O.bodies[0].state.vel[1];p1=O.bodies[0].state.pos[1]
		#print v1,p1
		O.run(10000,True)
		v2=O.bodies[0].state.vel[1];p2=O.bodies[0].state.pos[1]
		#print v2,p2
		self.assertTrue(abs(tv1-v1)<abs(tolerance*tv1) and abs(tp1-p1)<abs(tolerance*tp1))
		self.assertTrue(abs(tp2-p2)<abs(tolerance*tp2))
		#self.assertTrue(abs(tv2-v2)<abs(tolerance*tv2) and abs(tp2-p2)<abs(tolerance*tp2)) #velocity comparison disabled, see comment above
