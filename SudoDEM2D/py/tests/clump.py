
'''
Various computations affected by the periodic boundary conditions.
'''

import unittest
import random
from sudodem.wrapper import *
from sudodem._customConverters import *
from sudodem import utils
from sudodem import *
from math import *

try:
	from minieigen import *
except ImportError:
	from miniEigen import *

class TestSimpleClump(unittest.TestCase):
	"Test things on a simple clump composed of 2 disks."
	def setUp(self):
		O.reset()
		r1,r2,p0,p1=1,.5,Vector3.Zero,Vector3(0,0,3)
		self.idC,(self.id1,self.id2)=O.bodies.appendClumped([
			utils.disk(p0,r1),
			utils.disk(p1,r2)
		])
	def testConsistency(self):
		"Clump: ids and flags consistency"
		b1,b2,bC=[O.bodies[id] for id in (self.id1,self.id2,self.idC)]
		self.assertEqual(b1.clumpId,bC.id)
		self.assertEqual(b2.clumpId,bC.id)
		self.assertEqual(bC.clumpId,bC.id)
		self.assert_(bC.isClump)
		self.assert_(b1.isClumpMember)
		self.assert_(b2.isClumpMember)
		self.assert_(not bC.bounded)
	def testStaticProperties(self):
		"Clump: mass, centroid, intertia"
		b1,b2,bC=[O.bodies[id] for id in (self.id1,self.id2,self.idC)]
		# mass
		self.assertEqual(bC.state.mass,b1.state.mass+b2.state.mass)
		# centroid
		S=b1.state.mass*b1.state.pos+b2.state.mass*b2.state.pos
		c=S/bC.state.mass
		self.assertAlmostEqual(bC.state.pos[0],c[0]);
		self.assertAlmostEqual(bC.state.pos[1],c[1]);
		self.assertAlmostEqual(bC.state.pos[2],c[2]);
		# inertia
		i1,i2=(8./15)*pi*b1.material.density*b1.shape.radius**5, (8./15)*pi*b2.material.density*b2.shape.radius**5 # inertia of disks
		iMax=i1+i2+b1.state.mass*(b1.state.pos-c).norm()**2+b2.state.mass*(b2.state.pos-c).norm()**2 # minimum principal inertia
		iMin=i1+i2 # perpendicular to the
		# the order of bC.state.inertia is arbitrary (though must match the orientation)
		iC=list(bC.state.inertia); iC.sort()
		self.assertAlmostEqual(iC[0],iMin)
		self.assertAlmostEqual(iC[1],iMax)
		self.assertAlmostEqual(iC[2],iMax)
		# check orientation...?
		#self.assertAlmostEqual
	def testVelocity(self):
		"Clump: velocities of member assigned by NewtonIntegrator"
		s1,s2,sC=[O.bodies[id].state for id in (self.id1,self.id2,self.idC)]
		O.dt=0
		sC.vel=(1.,.2,.4)
		sC.angVel=(0,.4,.1)
		O.engines=[NewtonIntegrator()]; O.step() # update velocities
		# linear velocities
		self.assertEqual(s1.vel,sC.vel+sC.angVel.cross(s1.pos-sC.pos))
		self.assertEqual(s2.vel,sC.vel+sC.angVel.cross(s2.pos-sC.pos))
		# angular velocities
		self.assertEqual(s1.angVel,sC.angVel);
		self.assertEqual(s2.angVel,sC.angVel);

