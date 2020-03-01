# encoding: utf-8
# 2009 © Václav Šmilauer <eudoxos@arcig.cz>

"""
This test module covers python/c++ transitions, for both classes deriving from Serializable,
but also for other classes that we wrap (like miniEigen).
"""

import unittest
from sudodem.wrapper import *
from sudodem._customConverters import *
from math import *
from sudodem import system
from sudodem import *

try:
	from minieigen import *
except ImportError:
	from miniEigen import *

allClasses=system.childClasses('Serializable')

class TestObjectInstantiation(unittest.TestCase):
	def setUp(self):
		pass # no setup needed for tests here
	def testClassCtors(self):
		"Core: correct types are instantiated"
		# correct instances created with Foo() syntax
		for r in allClasses:
			obj=eval(r)();
			self.assert_(obj.__class__.__name__==r,'Failed for '+r)
	def testRootDerivedCtors_attrs_few(self):
		"Core: class ctor's attributes"
		# attributes passed when using the Foo(attr1=value1,attr2=value2) syntax
		gm=Shape(wire=True); self.assert_(gm.wire==True)
	def testDispatcherCtor(self):
		"Core: dispatcher ctors with functors"
		# dispatchers take list of their functors in the ctor
		# same functors are collapsed in one
		cld1=LawDispatcher([Law2_ScGeom_FrictPhys_CundallStrack(),Law2_ScGeom_FrictPhys_CundallStrack()]); self.assert_(len(cld1.functors)==1)
		# two different make two different, right?
		cld2=LawDispatcher([Law2_ScGeom_FrictPhys_CundallStrack(),Law2_ScGeom6D_CohFrictPhys_CohesionMoment()]); self.assert_(len(cld2.functors)==2)
	def testInteractionLoopCtor(self):
		"Core: InteractionLoop special ctor"
		# InteractionLoop takes 3 lists
		id=InteractionLoop([Ig2_Facet_Disk_ScGeom(),Ig2_Disk_Disk_ScGeom()],[Ip2_FrictMat_FrictMat_FrictPhys()],[Law2_ScGeom_FrictPhys_CundallStrack()],)
		self.assert_(len(id.geomDispatcher.functors)==2)
		self.assert_(id.geomDispatcher.__class__==IGeomDispatcher().__class__)
		self.assert_(id.physDispatcher.functors[0].__class__==Ip2_FrictMat_FrictMat_FrictPhys().__class__)
		self.assert_(id.lawDispatcher.functors[0].__class__==Law2_ScGeom_FrictPhys_CundallStrack().__class__)
	def testParallelEngineCtor(self):
		"Core: ParallelEngine special ctor"
		pe=ParallelEngine([InsertionSortCollider(),[BoundDispatcher(),ForceResetter()]])
		self.assert_(pe.slaves[0].__class__==InsertionSortCollider().__class__)
		self.assert_(len(pe.slaves[1])==2)
		pe.slaves=[]
		self.assert_(len(pe.slaves)==0)
	##
	## testing incorrect operations that should raise exceptions
	##
	def testWrongFunctorType(self):
		"Core: dispatcher and functor type mismatch is detected"
		# dispatchers accept only correct functors
		self.assertRaises(TypeError,lambda: LawDispatcher([Bo1_Disk_Aabb()]))
	def testInvalidAttr(self):
		'Core: invalid attribute access raises AttributeError'
		# accessing invalid attributes raises AttributeError
		self.assertRaises(AttributeError,lambda: Disk(attributeThatDoesntExist=42))
		self.assertRaises(AttributeError,lambda: Disk().attributeThatDoesntExist)
	##
	## attribute flags
	##
	def testTriggerPostLoad(self):
		'Core: Attr::triggerPostLoad'
		# TranslationEngine normalizes translationAxis automatically
		# anything else could be tested
		te=TranslationEngine();
		te.translationAxis=(0,2,0)
		self.assert_(te.translationAxis==(0,1,0))
	def testHidden(self):
		'Core: Attr::hidden'
		# hidden attributes are not wrapped in python at all
		self.assert_(not hasattr(Interaction(),'iterLastSeen'))
	def testNoSave(self):
		'Core: Attr::noSave'
		# update bound of the particle
		O.bodies.append(utils.disk((0,0,0),1))
		O.engines=[InsertionSortCollider([Bo1_Disk_Aabb()]),NewtonIntegrator()]
		O.step()
		O.saveTmp(quiet=True)
		mn0=Vector3(O.bodies[0].bound.min)
		O.reload(quiet=True)
		mn1=Vector3(O.bodies[0].bound.min)
		# check that the minimum is not saved
		self.assert_(not isnan(mn0[0]))
		self.assert_(isnan(mn1[0]))
	def _testReadonly(self):
		'Core: Attr::readonly'
		self.assertRaises(AttributeError,lambda: setattr(Body(),'id',3))

class TestEigenWrapper(unittest.TestCase):
	def assertSeqAlmostEqual(self,v1,v2):
		"floating-point comparison of vectors/quaterions"
		self.assertEqual(len(v1),len(v2));
		for i in range(len(v1)): self.assertAlmostEqual(v1[i],v2[i],msg='Component '+str(i)+' of '+str(v1)+' and '+str(v2))
	def testVector2(self):
		"Math: Vector2 operations"
		v=Vector2(1,2); v2=Vector2(3,4)
		self.assert_(v+v2==Vector2(4,6))
		self.assert_(Vector2().UnitX.dot(Vector2().UnitY)==0)
		self.assert_(Vector2().Zero.norm()==0)
	def testVector3(self):
		"Math: Vector3 operations"
		v=Vector3(3,4,5); v2=Vector3(3,4,5)
		self.assert_(v[0]==3 and v[1]==4 and v[2]==5)
		self.assert_(v.squaredNorm()==50)
		self.assert_(v==(3,4,5)) # comparison with list/tuple
		self.assert_(v==[3,4,5])
		self.assert_(v==v2)
		x,y,z,one=Vector3().UnitX,Vector3().UnitY,Vector3().UnitZ,Vector3().Ones
		self.assert_(x+y+z==one)
		self.assert_(x.dot(y)==0)
		self.assert_(x.cross(y)==z)
	def testQuaternion(self):
		"Math: Quaternion operations"
		# construction
		q1=Quaternion((0,0,1),pi/2)
		q2=Quaternion(Vector3(0,0,1),pi/2)
		q1==q2
		x,y,z,one=Vector3().UnitX,Vector3().UnitY,Vector3().UnitZ,Vector3().Ones
		self.assertSeqAlmostEqual(q1*x,y)
		self.assertSeqAlmostEqual(q1*q1*x,-x)
		self.assertSeqAlmostEqual(q1*q1.conjugate(),Quaternion().Identity)
		self.assertSeqAlmostEqual(q1.toAxisAngle()[0],(0,0,1))
		self.assertAlmostEqual(q1.toAxisAngle()[1],pi/2)
	def testMatrix3(self):
		"Math: Matrix3 operations"
		#construction
		m1=Matrix3(1,0,0,0,1,0,0,0,1)
		# comparison
		self.assert_(m1==Matrix3().Identity)
		# rotation matrix from quaternion
		m1=Quaternion(Vector3(0,0,1),pi/2).toRotationMatrix()
		# multiplication with vectors
		self.assertSeqAlmostEqual(m1*Vector3().UnitX,Vector3().UnitY)
		# determinant
		m2=Matrix3(-2,2,-3,-1,1,3,2,0,-1)
		self.assertEqual(m2.determinant(),18)
		# inverse
		inv=Matrix3(-0.055555555555556,0.111111111111111,0.5,0.277777777777778,0.444444444444444,0.5,-0.111111111111111,0.222222222222222,0.0)
		m2inv=m2.inverse()
		self.assertSeqAlmostEqual(m2inv,inv)
		# matrix-matrix multiplication
		self.assertSeqAlmostEqual(Matrix3().Identity*Matrix3().Identity,Matrix3().Identity)
		m3=Matrix3(1,2,3,4,5,6,-1,0,3)
		m33=m3*m3
		self.assertSeqAlmostEqual(m33,Matrix3(6,12,24,18,33,60,-4,-2,6))

	# not really wm3 thing, but closely related
	# no way to test this currently, as State::se3 is not serialized (State::pos and State::ori are serialized instead...)
	#def testSe3Conversion(self):
	#	return
	#	pp=State()
	#	pp.se3=(Vector3().Zero,Quaternion().Identity)
	#	self.assert_(pp['se3'][0]==Vector3().Zero)
	#	self.assert_(pp['se3'][1]==Quaternion().Identity)
	#	pp.se3=((1,2,3),Quaternion((1,1,1),pi/4))
	#	self.assert_(pp['se3'][0]==(1,2,3))
	#	self.assert_(pp['se3'][0]==pp.pos)
	#	self.assert_(pp['se3'][1]==Quaternion((1,1,1),pi/4))
	#	self.assert_(pp['se3'][1]==pp.ori)


