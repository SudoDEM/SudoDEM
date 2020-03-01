# encoding: utf-8
# 2010 © Václav Šmilauer <eudoxos@arcig.cz>

'''
Various computations affected by the periodic boundary conditions.
'''

import unittest
import random,math
from sudodem.wrapper import *
from sudodem._customConverters import *
from sudodem import utils
from sudodem import *

try:
	from minieigen import *
except ImportError:
	from miniEigen import *

class TestPBC(unittest.TestCase):
	# prefix test names with PBC:
	def setUp(self):
		O.reset(); O.periodic=True;
		O.cell.setBox(2.5,2.5,3)
		self.cellDist=Vector3i(0,0,10) # how many cells away we go
		self.relDist=Vector3(0,.999999999999999999,0) # rel position of the 2nd ball within the cell
		self.initVel=Vector3(0,0,5)
		O.bodies.append(utils.disk((1,1,1),.5))
		self.initPos=Vector3([O.bodies[0].state.pos[i]+self.relDist[i]+self.cellDist[i]*O.cell.refSize[i] for i in (0,1,2)])
		O.bodies.append(utils.disk(self.initPos,.5))
		#print O.bodies[1].state.pos
		O.bodies[1].state.vel=self.initVel
		O.engines=[NewtonIntegrator(warnNoForceReset=False)]
		O.cell.velGrad=Matrix3(0,0,0, 0,0,0, 0,0,-1)
		O.dt=0 # do not change positions with dt=0 in NewtonIntegrator, but still update velocities from velGrad
	def testVelGrad(self):
		'PBC: velGrad changes hSize, accumulates in trsf'
		O.dt=1e-3
		hSize,trsf=O.cell.hSize,Matrix3.Identity
		hSize0=hSize
		O.step()
		for i in range(0,10):
			O.step(); hSize+=O.dt*O.cell.velGrad*hSize; trsf+=O.dt*O.cell.velGrad*trsf
		for i in range(0,len(O.cell.hSize)):
			self.assertAlmostEqual(hSize[i],O.cell.hSize[i])
			self.assertAlmostEqual(trsf[i],O.cell.trsf[i])
	def testTrsfChange(self):
		'PBC: chaing trsf changes hSize0, but does not modify hSize'
		O.dt=1e-2
		O.step()
		O.cell.trsf=Matrix3.Identity
		for i in range(0,len(O.cell.hSize)):
			self.assertAlmostEqual(O.cell.hSize0[i],O.cell.hSize[i])
	def testDegenerate(self):
		"PBC: degenerate cell raises exception"
		self.assertRaises(RuntimeError,lambda: setattr(O.cell,'hSize',Matrix3(1,0,0, 0,0,0, 0,0,1)))
	def testSetBox(self):
		"PBC: setBox modifies hSize correctly"
		O.cell.setBox(2.55,11,45)
		self.assert_(O.cell.hSize==Matrix3(2.55,0,0, 0,11,0, 0,0,45));
	def testHomotheticResizeVel(self):
		"PBC: homothetic cell deformation adjusts particle velocity "
		O.dt=1e-5
		O.step()
		s1=O.bodies[1].state
		self.assertAlmostEqual(s1.vel[2],self.initVel[2]+self.initPos[2]*O.cell.velGrad[2,2])
	def testScGeomIncidentVelocity(self):
		"PBC: ScGeom computes incident velocity correctly"
		O.run(2,1)
		O.engines=[InteractionLoop([Ig2_Disk_Disk_ScGeom()],[Ip2_FrictMat_FrictMat_FrictPhys()],[])]
		i=utils.createInteraction(0,1)
		self.assertEqual(self.initVel,i.geom.incidentVel(i,avoidGranularRatcheting=True))
		self.assertEqual(self.initVel,i.geom.incidentVel(i,avoidGranularRatcheting=False))
		self.assertAlmostEqual(self.relDist[1],1-i.geom.penetrationDepth)
	def testL3GeomIncidentVelocity(self):
		"PBC: L3Geom computes incident velocity correctly"
		O.step()
		O.engines=[ForceResetter(),InteractionLoop([Ig2_Disk_Disk_L3Geom()],[Ip2_FrictMat_FrictMat_FrictPhys()],[Law2_L3Geom_FrictPhys_ElPerfPl(noBreak=True)]),NewtonIntegrator()]
		i=utils.createInteraction(0,1)
		O.dt=1e-10; O.step() # tiny timestep, to not move the normal too much
		self.assertAlmostEqual(self.initVel.norm(),(i.geom.u/O.dt).norm())
	def testKineticEnergy(self):
		"PBC: utils.kineticEnergy considers only fluctuation velocity, not the velocity gradient"
		O.step() # updates velocity with homotheticCellResize
		# ½(mv²+ωIω)
		# #0 is still, no need to add it; #1 has zero angular velocity
		# we must take self.initVel since O.bodies[1].state.vel now contains the homothetic resize which utils.kineticEnergy is supposed to compensate back
		Ek=.5*O.bodies[1].state.mass*self.initVel.squaredNorm()
		self.assertAlmostEqual(Ek,utils.kineticEnergy())



