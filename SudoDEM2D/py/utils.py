# encoding: utf-8
#
# utility functions for sudodem
#
# 2008-2009 © Václav Šmilauer <eudoxos@arcig.cz>

"""Heap of functions that don't (yet) fit anywhere else.

Devs: please DO NOT ADD more functions here, it is getting too crowded!
"""

import math,random,doctest,numpy
from sudodem import geom
from sudodem import *
from sudodem.wrapper import *
try: # use psyco if available
	import psyco
	psyco.full()
except ImportError: pass

try:
	from minieigen import *
except ImportError:
	from miniEigen import *


# c++ implementations for performance reasons
from sudodem._utils import *

def saveVars(mark='',loadNow=True,**kw):
	"""Save passed variables into the simulation so that it can be recovered when the simulation is loaded again.

	For example, variables *a*, *b* and *c* are defined. To save them, use::

		>>> saveVars('something',a=1,b=2,c=3)
		>>> from sudodem.params.something import *
		>>> a,b,c
		(1, 2, 3)

	those variables will be save in the .xml file, when the simulation itself is saved. To recover those variables once the .xml is loaded again, use
	``loadVars('something')``and they will be defined in the sudodem.params.\ *mark* module. The *loadNow* parameter calls :yref:`sudodem.utils.loadVars`
	after saving automatically. If 'something' already exists, given variables will be inserted.
	"""
	import cPickle
	try:
		d=cPickle.loads(Omega().tags['pickledPythonVariablesDictionary'+mark])	#load dictionary d
		for key in kw.keys():
			d[key]=kw[key]							#insert new variables into d
	except KeyError:
		d = kw
	Omega().tags['pickledPythonVariablesDictionary'+mark]=cPickle.dumps(d)
	if loadNow: loadVars(mark)


def loadVars(mark=None):
	"""Load variables from :yref:`sudodem.utils.saveVars`, which are saved inside the simulation.
	If ``mark==None``, all save variables are loaded. Otherwise only those with
	the mark passed."""
	import cPickle, types, sys, warnings
	def loadOne(d,mark=None):
		"""Load given dictionary into a synthesized module sudodem.params.name (or sudodem.params if *name* is not given). Update sudodem.params.__all__ as well."""
		import sudodem.params
		if mark:
			if mark in sudodem.params.__dict__: warnings.warn('Overwriting sudodem.params.%s which already exists.'%mark)
			modName='sudodem.params.'+mark
			mod=types.ModuleType(modName)
			mod.__dict__.update(d)
			mod.__all__=list(d.keys()) # otherwise params starting with underscore would not be imported
			sys.modules[modName]=mod
			sudodem.params.__all__.append(mark)
			sudodem.params.__dict__[mark]=mod
		else:
			sudodem.params.__all__+=list(d.keys())
			sudodem.params.__dict__.update(d)
	if mark!=None:
		d=cPickle.loads(Omega().tags['pickledPythonVariablesDictionary'+mark])
		loadOne(d,mark)
	else: # load everything one by one
		for m in Omega().tags.keys():
			if m.startswith('pickledPythonVariablesDictionary'):
				loadVars(m[len('pickledPythonVariableDictionary')+1:])


def DiskPWaveTimeStep(radius,density,young):
	r"""Compute P-wave critical timestep for a single (presumably representative) disk, using formula for P-Wave propagation speed $\Delta t_{c}=\frac{r}{\sqrt{E/\rho}}$.
	If you want to compute minimum critical timestep for all disks in the simulation, use :yref:`sudodem.utils.PWaveTimeStep` instead.

	>>> DiskPWaveTimeStep(1e-3,2400,30e9)
	2.8284271247461903e-07
	"""
	from math import sqrt
	return radius/sqrt(young/density)

def randomColor():
	"""Return random Vector3 with each component in interval 0…1 (uniform distribution)"""
	return Vector3(random.random(),random.random(),random.random())

def typedEngine(name):
	"""Return first engine from current O.engines, identified by its type (as string). For example:

	>>> from sudodem import utils
	>>> O.engines=[InsertionSortCollider(),NewtonIntegrator(),GravityEngine()]
	>>> utils.typedEngine("NewtonIntegrator") == O.engines[1]
	True
	"""
	return [e for e in Omega().engines if e.__class__.__name__==name][0]

def defaultMaterial():
	"""Return default material, when creating bodies with :yref:`sudodem.utils.disk` and friends, material is unspecified and there is no shared material defined yet. By default, this function returns::

		.. code-block:: python

		  FrictMat(density=1e3,young=1e7,poisson=.3,frictionAngle=.5,label='defaultMat')
	"""
	return FrictMat(density=1e3,young=1e7,poisson=.3,frictionAngle=.5,label='defaultMat')

def _commonBodySetup(b,volume,geomInertia,material,pos,noBound=False,resetState=True,dynamic=None,fixed=False):
	"""Assign common body parameters."""
	if isinstance(material,int):
		if material<0 and len(O.materials)==0: O.materials.append(defaultMaterial());
		b.mat=O.materials[material]
	elif isinstance(material,str): b.mat=O.materials[material]
	elif isinstance(material,Material): b.mat=material
	elif callable(material): b.mat=material()
	else: raise TypeError("The 'material' argument must be None (for defaultMaterial), string (for shared material label), int (for shared material id) or Material instance.");
	## resets state (!!)
	if resetState: b.state=b.mat.newAssocState()
	mass=volume*b.mat.density
	b.state.mass,b.state.inertia=mass,geomInertia*b.mat.density
	b.state.pos=b.state.refPos=pos
	b.bounded=(not noBound)
	if dynamic!=None:
		import warnings
		warnings.warn('dynamic=%s is deprecated, use fixed=%s instead'%(str(dynamic),str(not dynamic)),category=DeprecationWarning,stacklevel=2)
		fixed=not dynamic
	b.state.blockedDOFs=('xyZ' if fixed else '')


def disk(center,radius,z_dim=1,dynamic=None,fixed=False,wire=False,color=None,highlight=False,material=-1,mask=1):
	"""Create disk with given parameters; mass and inertia computed automatically.

	Last assigned material is used by default (*material* = -1), and utils.defaultMaterial() will be used if no material is defined at all.

	:param Vector2 center: center
	:param float radius: radius
	:param float dynamic: deprecated, see "fixed"
	:param float fixed: generate the body with all DOFs blocked?
	:param material:
		specify :yref:`Body.material`; different types are accepted:
			* int: O.materials[material] will be used; as a special case, if material==-1 and there is no shared materials defined, utils.defaultMaterial() will be assigned to O.materials[0]
			* string: label of an existing material that will be used
			* :yref:`Material` instance: this instance will be used
			* callable: will be called without arguments; returned Material value will be used (Material factory object, if you like)
	:param int mask: :yref:`Body.mask` for the body
	:param wire: display as wire disk?
	:param highlight: highlight this body in the viewer?
	:param Vector2-or-None: body's color, as normalized RGB; random color will be assigned if ``None``.

	:return:
		A Body instance with desired characteristics.


	Creating default shared material if none exists neither is given::

		>>> O.reset()
		>>> from sudodem import utils
		>>> len(O.materials)
		0
		>>> s0=utils.disk([2,0,0],1)
		>>> len(O.materials)
		1

	Instance of material can be given::

		>>> s1=utils.disk([0,0,0],1,wire=False,color=(0,1,0),material=ElastMat(young=30e9,density=2e3))
		>>> s1.shape.wire
		False
		>>> s1.shape.color
		Vector2(0,1,0)
		>>> s1.mat.density
		2000.0

	Material can be given by label::

		>>> O.materials.append(FrictMat(young=10e9,poisson=.11,label='myMaterial'))
		1
		>>> s2=utils.disk([0,0,2],1,material='myMaterial')
		>>> s2.mat.label
		'myMaterial'
		>>> s2.mat.poisson
		0.11

	Finally, material can be a callable object (taking no arguments), which returns a Material instance.
	Use this if you don't call this function directly (for instance, through sudodem.pack.randomDensePack), passing
	only 1 *material* parameter, but you don't want material to be shared.

	For instance, randomized material properties can be created like this:

		>>> import random
		>>> def matFactory(): return ElastMat(young=1e10*random.random(),density=1e3+1e3*random.random())
		...
		>>> s3=utils.disk([0,2,0],1,material=matFactory)
		>>> s4=utils.disk([1,2,0],1,material=matFactory)

	"""
	b=Body()
	b.shape=Disk(radius=radius,ref_radius=radius,color=color if color else randomColor(),wire=wire,highlight=highlight)
	V= math.pi*radius**2*z_dim
	geomInert=0.5*V*radius**2
	_commonBodySetup(b,V,geomInert,material,pos=center,dynamic=dynamic,fixed=fixed)
	b.aspherical=False
	b.mask=mask
	return b

def box(center,extents,orientation=Quaternion(1,0,0,0),dynamic=None,fixed=False,wire=False,color=None,highlight=False,material=-1,mask=1):
	"""Create box (cuboid) with given parameters.

	:param Vector2 extents: half-sizes along x,y,z axes

	See :yref:`sudodem.utils.disk`'s documentation for meaning of other parameters."""
	b=Body()
	b.shape=Box(extents=extents,color=color if color else randomColor(),wire=wire,highlight=highlight)
	V=8*extents[0]*extents[1]*extents[2]
	geomInert=Vector3(4*(extents[1]**2+extents[2]**2),4*(extents[0]**2+extents[2]**2),4*(extents[0]**2+extents[1]**2))
	_commonBodySetup(b,V,geomInert,material,pos=center,dynamic=dynamic,fixed=fixed)
	b.state.ori=orientation
	b.mask=mask
	b.aspherical=True
	return b


def chainedCylinder(begin=Vector3(0,0,0),end=Vector3(1.,0.,0.),radius=0.2,dynamic=None,fixed=False,wire=False,color=None,highlight=False,material=-1,mask=1):
	"""Create and connect a chainedCylinder with given parameters. The shape generated by repeted calls of this function is the Minkowski sum of polyline and disk.

	:param Real radius: radius of disk in the Minkowski sum.
	:param Vector3 begin: first point positioning the line in the Minkowski sum
	:param Vector3 last: last point positioning the line in the Minkowski sum

	In order to build a correct chain, last point of element of rank N must correspond to first point of element of rank N+1 in the same chain (with some tolerance, since bounding boxes will be used to create connections.

	:return: Body object with the :yref:`ChainedCylinder` :yref:`shape<Body.shape>`.
	"""
	segment=end-begin
	b=Body()
	b.shape=ChainedCylinder(radius=radius,length=segment.norm(),color=color if color else randomColor(),wire=wire,highlight=highlight)
	b.shape.segment=segment
	V=2*(4./3)*math.pi*radius**3
	geomInert=(2./5.)*V*radius**2+b.shape.length*b.shape.length*2*(4./3)*math.pi*radius**3
	b.state=ChainedState()
	b.state.addToChain(O.bodies.append(b))
	_commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=begin,resetState=False,dynamic=dynamic,fixed=fixed)
	b.mask=mask
	b.bound=Aabb(color=[0,1,0])
	b.state.ori.setFromTwoVectors(Vector3(0.,0.,1.),segment)
	if (end == begin): b.state.ori = Quaternion((1,0,0),0)
	return b

def gridNode(center,radius,dynamic=None,fixed=False,wire=False,color=None,highlight=False,material=-1):
	b=Body()
	b.shape=GridNode(radius=radius,color=color if color else randomColor(),wire=wire,highlight=highlight)
	#V=(4./3)*math.pi*radius**3	# will be overriden by the connection
	V=0.
	geomInert=(2./5.)*V*radius**2	# will be overriden by the connection
	_commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=center,dynamic=dynamic,fixed=fixed)
	b.aspherical=False
	b.bounded=False
	b.mask=0	#avoid contact detection with the nodes. Manual interaction will be set for them in "gridConnection" below.
	return b

def gridConnection(id1,id2,radius,wire=False,color=None,highlight=False,material=-1,mask=1,cellDist=None):
	b=Body()
	b.shape=GridConnection(radius=radius,color=color if color else randomColor(),wire=wire,highlight=highlight)
	sph1=O.bodies[id1] ; sph2=O.bodies[id2]
	i=createInteraction(id1,id2)
	nodeMat=sph1.material
	b.shape.node1=sph1 ; b.shape.node2=sph2
	sph1.shape.addConnection(b) ; sph2.shape.addConnection(b)
	if(O.periodic):
		if(cellDist!=None):
			i.cellDist=cellDist
		segt=sph2.state.pos + O.cell.hSize*i.cellDist - sph1.state.pos
	else: segt=sph2.state.pos - sph1.state.pos
	L=segt.norm()
	V=0.5*L*math.pi*radius**2
	geomInert=(2./5.)*V*radius**2
	_commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=sph1.state.pos,dynamic=False,fixed=True)
	sph1.state.mass = sph1.state.mass + V*nodeMat.density
	sph2.state.mass = sph2.state.mass + V*nodeMat.density
	for k in [0,1,2]:
		sph1.state.inertia[k] = sph1.state.inertia[k] + geomInert*nodeMat.density
		sph2.state.inertia[k] = sph2.state.inertia[k] + geomInert*nodeMat.density
	b.aspherical=False
	if O.periodic:
		i.phys.unp= -(sph2.state.pos + O.cell.hSize*i.cellDist - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius
		b.shape.periodic=True
		b.shape.cellDist=i.cellDist
	else:
		i.phys.unp= -(sph2.state.pos - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius
	i.geom.connectionBody=b
	I=math.pi*(2.*radius)**4/64.
	E=nodeMat.young
	i.phys.kn=E*math.pi*(radius**2)/L
	i.phys.kr=E*I/L
	i.phys.ks=12.*E*I/(L**3)
	G=E/(2.*(1+nodeMat.poisson))
	i.phys.ktw=2.*I*G/L
	b.mask=mask
	return b


def wall(position,axis,sense=0,color=None,material=-1,mask=1):
	"""Return ready-made wall body.

	:param float-or-Vector3 position: center of the wall. If float, it is the position along given axis, the other 2 components being zero
	:param ∈{0,1,2} axis: orientation of the wall normal (0,1,2) for x,y,z (sc. planes yz, xz, xy)
	:param ∈{-1,0,1} sense: sense in which to interact (0: both, -1: negative, +1: positive; see :yref:`Wall`)

	See :yref:`sudodem.utils.disk`'s documentation for meaning of other parameters."""
	b=Body()
	b.shape=Wall(sense=sense,axis=axis,color=color if color else randomColor())
	if isinstance(position,(int,long,float)):
		pos2=Vector2(0,0); pos2[axis]=position
	else: pos2=position
	_commonBodySetup(b,0,0,material,pos=pos2,fixed=True)
	b.aspherical=False # wall never moves dynamically
	b.mask=mask
	return b

def fwall(vertex1,vertex2,dynamic=None,fixed=True,color=None,highlight=False,noBound=False,material=-1,mask=1,chain=-1):
	"""Create fwall with given parameters.

	:param Vector2 vertex1: coordinates of vertex1 in the global coordinate system.
	:param Vector2 vertex2: coordinates of vertex2 in the global coordinate system.
	:param bool noBound: set :yref:`Body.bounded`
	:param Vector3-or-None color: color of the facet; random color will be assigned if ``None``.

	See :yref:`sudodem.utils.disk`'s documentation for meaning of other parameters."""
	b=Body()
	center=0.5*(Vector2(vertex1) + Vector2(vertex2))
	v1 = Vector2(vertex1) - center
	v2 = Vector2(vertex2) - center
	b.shape=Fwall(color=color if color else randomColor(),highlight=highlight,vertex1 = v1, vertex2 = v2)
	_commonBodySetup(b,0,0,material,noBound=noBound,pos=center,fixed=fixed)
	b.aspherical=False # mass and inertia are 0 anyway; fell free to change to ``True`` if needed
	b.mask=mask
	b.chain=chain
	return b

def facet(vertices,dynamic=None,fixed=True,wire=True,color=None,highlight=False,noBound=False,material=-1,mask=1,chain=-1):
	"""Create facet with given parameters.

	:param [Vector3,Vector3,Vector3] vertices: coordinates of vertices in the global coordinate system.
	:param bool wire: if ``True``, facets are shown as skeleton; otherwise facets are filled
	:param bool noBound: set :yref:`Body.bounded`
	:param Vector3-or-None color: color of the facet; random color will be assigned if ``None``.

	See :yref:`sudodem.utils.disk`'s documentation for meaning of other parameters."""
	b=Body()
	center=inscribedCircleCenter(vertices[0],vertices[1],vertices[2])
	vertices=Vector3(vertices[0])-center,Vector3(vertices[1])-center,Vector3(vertices[2])-center
	b.shape=Facet(color=color if color else randomColor(),wire=wire,highlight=highlight,vertices=vertices)
	_commonBodySetup(b,0,Vector3(0,0,0),material,noBound=noBound,pos=center,fixed=fixed)
	b.aspherical=False # mass and inertia are 0 anyway; fell free to change to ``True`` if needed
	b.mask=mask
	b.chain=chain
	return b

def tetraPoly(vertices,dynamic=True,fixed=False,wire=True,color=None,highlight=False,noBound=False,material=-1,mask=1,chain=-1):
	"""Create tetrahedron (actually simple Polyhedra) with given parameters.

	:param [Vector3,Vector3,Vector3,Vector3] vertices: coordinates of vertices in the global coordinate system.

	See :yref:`sudodem.utils.disk`'s documentation for meaning of other parameters."""
	b=Body()
	b.shape = Polyhedra(v=vertices,color=color if color else randomColor(),wire=wire,highlight=highlight)
	volume = b.shape.GetVolume()
	inertia = b.shape.GetInertia()
	center = b.shape.GetCentroid()
	_commonBodySetup(b,volume,inertia,material,noBound=noBound,pos=center,fixed=fixed)
	b.aspherical=False
	b.state.ori = b.shape.GetOri()
	b.mask=mask
	b.chain=chain
	return b

def tetra(vertices,strictCheck=True,dynamic=True,fixed=False,wire=True,color=None,highlight=False,noBound=False,material=-1,mask=1,chain=-1):
	"""Create tetrahedron with given parameters.

	:param [Vector3,Vector3,Vector3,Vector3] vertices: coordinates of vertices in the global coordinate system.
	:param bool strictCheck: checks vertices order, raise RuntimeError for negative volume

	See :yref:`sudodem.utils.disk`'s documentation for meaning of other parameters."""
	b=Body()
	center = .25*sum(vertices,Vector3.Zero)
	volume = TetrahedronSignedVolume(vertices)
	if volume < 0:
		if strictCheck:
			raise (RuntimeError, "tetra: wrong order of vertices")
		temp = vertices[3]
		vertices[3] = vertices[2]
		vertices[2] = temp
		volume = TetrahedronSignedVolume(vertices)
	assert(volume>0)
	b.shape = Tetra(v=vertices,color=color if color else randomColor(),wire=wire,highlight=highlight)
	# modifies pos, ori and inertia
	ori = TetrahedronWithLocalAxesPrincipal(b)
	_commonBodySetup(b,volume,b.state.inertia,material,noBound=noBound,pos=center,fixed=fixed)
	b.state.ori = b.state.refOri = ori
	b.aspherical = True
	b.mask = mask
	b.chain = chain
	return b

def polyhedron(vertices,dynamic=True,fixed=False,wire=True,color=None,highlight=False,noBound=False,material=-1,mask=1,chain=-1):
	"""Create polyhedron with given parameters.

	:param [Vector3,Vector3,Vector3,Vector3] vertices: coordinates of vertices in the global coordinate system.

	See :yref:`sudodem.utils.disk`'s documentation for meaning of other parameters."""
	b=Body()
	b.shape = Polyhedra(v=vertices,color=color if color else randomColor(),wire=wire,highlight=highlight)
	volume = b.shape.GetVolume()
	inertia = b.shape.GetInertia()
	center = b.shape.GetCentroid()
	_commonBodySetup(b,volume,inertia,material,noBound=noBound,pos=center,fixed=fixed)
	b.aspherical=False
	b.state.ori = b.shape.GetOri()
	b.mask=mask
	b.chain=chain
	return b






#def setNewVerticesOfFacet(b,vertices):
#	center = inscribedCircleCenter(vertices[0],vertices[1],vertices[2])
#	vertices = Vector3(vertices[0])-center,Vector3(vertices[1])-center,Vector3(vertices[2])-center
#	b.shape.vertices = vertices
#	b.state.pos = center


def facetBox(*args,**kw):
	"|ydeprecated|"
	_deprecatedUtilsFunction('facetBox','geom.facetBox')
	return geom.facetBox(*args,**kw)

def facetCylinder(*args,**kw):
	"|ydeprecated|"
	_deprecatedUtilsFunction('facetCylinder','geom.facetCylinder')
	return geom.facetCylinder(*args,**kw)

def aabbWalls(extrema=None,thickness=0,oversizeFactor=1.5,**kw):
	"""Return 6 boxes that will wrap existing packing as walls from all sides;
	extrema are extremal points of the Aabb of the packing (will be calculated if not specified)
	thickness is wall thickness (will be 1/10 of the X-dimension if not specified)
	Walls will be enlarged in their plane by oversizeFactor.
	returns list of 6 wall Bodies enclosing the packing, in the order minX,maxX,minY,maxY,minZ,maxZ.
	"""
	walls=[]
	if not extrema: extrema=aabbExtrema()
	#if not thickness: thickness=(extrema[1][0]-extrema[0][0])/10.
	for axis in [0,1,2]:
		mi,ma=extrema
		center=[(mi[i]+ma[i])/2. for i in range(3)]
		extents=[.5*oversizeFactor*(ma[i]-mi[i]) for i in range(3)]
		extents[axis]=thickness/2.
		for j in [0,1]:
			center[axis]=extrema[j][axis]+(j-.5)*thickness
			walls.append(box(center=center,extents=extents,fixed=True,**kw))
			walls[-1].shape.wire=True
	return walls


def aabbDim(cutoff=0.,centers=False):
	"""Return dimensions of the axis-aligned bounding box, optionally with relative part *cutoff* cut away."""
	a=aabbExtrema(cutoff,centers)
	return (a[1][0]-a[0][0],a[1][1]-a[0][1],a[1][2]-a[0][2])

def aabbExtrema2d(pts):
	"""Return 2d bounding box for a sequence of 2-tuples."""
	inf=float('inf')
	min,max=[inf,inf],[-inf,-inf]
	for pt in pts:
		if pt[0]<min[0]: min[0]=pt[0]
		elif pt[0]>max[0]: max[0]=pt[0]
		if pt[1]<min[1]: min[1]=pt[1]
		elif pt[1]>max[1]: max[1]=pt[1]
	return tuple(min),tuple(max)

def perpendicularArea(axis):
	"""Return area perpendicular to given axis (0=x,1=y,2=z) generated by bodies
	for which the function consider returns True (defaults to returning True always)
	and which is of the type :yref:`Disk`.
	"""
	ext=aabbExtrema()
	other=((axis+1)%3,(axis+2)%3)
	return (ext[1][other[0]]-ext[0][other[0]])*(ext[1][other[1]]-ext[0][other[1]])

def fractionalBox(fraction=1.,minMax=None):
	"""retrurn (min,max) that is the original minMax box (or aabb of the whole simulation if not specified)
	linearly scaled around its center to the fraction factor"""
	if not minMax: minMax=aabbExtrema()
	half=[.5*(minMax[1][i]-minMax[0][i]) for i in [0,1,2]]
	return (tuple([minMax[0][i]+(1-fraction)*half[i] for i in [0,1,2]]),tuple([minMax[1][i]-(1-fraction)*half[i] for i in [0,1,2]]))


def randomizeColors(onlyDynamic=False):
	"""Assign random colors to :yref:`Shape::color`.

	If onlyDynamic is true, only dynamic bodies will have the color changed.
	"""
	for b in O.bodies:
		color=(random.random(),random.random(),random.random())
		if b.dynamic or not onlyDynamic: b.shape.color=color

def avgNumInteractions(cutoff=0.,skipFree=False,considerClumps=False):
	r"""Return average numer of interactions per particle, also known as *coordination number* $Z$. This number is defined as

	.. math:: Z=2C/N

	where $C$ is number of contacts and $N$ is number of particles. When clumps are present, number of particles is the sum of standalone disks plus the sum of clumps.
	Clumps are considered in the calculation if cutoff != 0 or skipFree = True. If cutoff=0 (default) and skipFree=False (default) one needs to set considerClumps=True to consider clumps in the calculation.

	With *skipFree*, particles not contributing to stable state of the packing are skipped, following equation (8) given in [Thornton2000]_:

	.. math:: Z_m=\frac{2C-N_1}{N-N_0-N_1}

	:param cutoff: cut some relative part of the sample's bounding box away.
	:param skipFree: see above.
	:param considerClumps: also consider clumps if cutoff=0 and skipFree=False; for further explanation see above.

"""
	if cutoff==0 and not skipFree and not considerClumps: return 2*O.interactions.countReal()*1./len(O.bodies)
	else:
		nums,counts=bodyNumInteractionsHistogram(aabbExtrema(cutoff))
		## CC is 2*C
		CC=sum([nums[i]*counts[i] for i in range(len(nums))]); N=sum(counts)
		if not skipFree: return CC*1./N if N>0 else float('nan')
		## find bins with 0 and 1 disks
		N0=0 if (0 not in nums) else counts[nums.index(0)]
		N1=0 if (1 not in nums) else counts[nums.index(1)]
		NN=N-N0-N1
		return (CC-N1)*1./NN if NN>0 else float('nan')

def plotNumInteractionsHistogram(cutoff=0.):
	"Plot histogram with number of interactions per body, optionally cutting away *cutoff* relative axis-aligned box from specimen margin."
	nums,counts=bodyNumInteractionsHistogram(aabbExtrema(cutoff))
	import pylab
	pylab.bar(nums,counts)
	pylab.title('Number of interactions histogram, average %g (cutoff=%g)'%(avgNumInteractions(cutoff),cutoff))
	pylab.xlabel('Number of interactions')
	pylab.ylabel('Body count')
	pylab.ion()
	pylab.show()

def plotDirections(aabb=(),mask=0,bins=20,numHist=True,noShow=False):
	"""Plot 3 histograms for distribution of interaction directions, in yz,xz and xy planes and
	(optional but default) histogram of number of interactions per body.

	:returns: If *noShow* is ``False``, displays the figure and returns nothing. If *noShow*, the figure object is returned without being displayed (works the same way as :yref:`sudodem.plot.plot`).
	"""
	import pylab,math
	from sudodem import utils
	for axis in [0,1,2]:
		d=utils.interactionAnglesHistogram(axis,mask=mask,bins=bins,aabb=aabb)
		fc=[0,0,0]; fc[axis]=1.
		subp=pylab.subplot(220+axis+1,polar=True);
		# 1.1 makes small gaps between values (but the column is a bit decentered)
		pylab.bar(d[0],d[1],width=math.pi/(1.1*bins),fc=fc,alpha=.7,label=['yz','xz','xy'][axis])
		#pylab.title(['yz','xz','xy'][axis]+' plane')
		pylab.text(.5,.25,['yz','xz','xy'][axis],horizontalalignment='center',verticalalignment='center',transform=subp.transAxes,fontsize='xx-large')
	if numHist:
		pylab.subplot(224,polar=False)
		nums,counts=utils.bodyNumInteractionsHistogram(aabb if len(aabb)>0 else utils.aabbExtrema())
		avg=sum([nums[i]*counts[i] for i in range(len(nums))])/(1.*sum(counts))
		pylab.bar(nums,counts,fc=[1,1,0],alpha=.7,align='center')
		pylab.xlabel('Interactions per body (avg. %g)'%avg)
		pylab.axvline(x=avg,linewidth=3,color='r')
		pylab.ylabel('Body count')
	if noShow: return pylab.gcf()
	else:
		pylab.ion()
		pylab.show()


def encodeVideoFromFrames(*args,**kw):
	"|ydeprecated|"
	_deprecatedUtilsFunction('utils.encodeVideoFromFrames','utils.makeVideo')
	return makeVideo(*args,**kw)

def makeVideo(frameSpec,out,renameNotOverwrite=True,fps=24,kbps=6000,bps=None):
	"""Create a video from external image files using `mencoder <http://www.mplayerhq.hu>`__. Two-pass encoding using the default mencoder codec (mpeg4) is performed, running multi-threaded with number of threads equal to number of OpenMP threads allocated for SudoDEM.

	:param frameSpec: wildcard | sequence of filenames. If list or tuple, filenames to be encoded in given order; otherwise wildcard understood by mencoder's mf:// URI option (shell wildcards such as ``/tmp/snap-*.png`` or and printf-style pattern like ``/tmp/snap-%05d.png``)
	:param str out: file to save video into
	:param bool renameNotOverwrite: if True, existing same-named video file will have -*number* appended; will be overwritten otherwise.
	:param int fps: Frames per second (``-mf fps=…``)
	:param int kbps: Bitrate (``-lavcopts vbitrate=…``) in kb/s
	"""
	import os,os.path,subprocess,warnings
	if bps!=None:
		warnings.warn('plot.makeVideo: bps is deprecated, use kbps instead (the significance is the same, but the name is more precise)',stacklevel=2,category=DeprecationWarning)
		kbps=bps
	if renameNotOverwrite and os.path.exists(out):
		i=0
		while(os.path.exists(out+"~%d"%i)): i+=1
		os.rename(out,out+"~%d"%i); print("Output file `%s' already existed, old file renamed to `%s'"%(out,out+"~%d"%i))
	if isinstance(frameSpec,list) or isinstance(frameSpec,tuple): frameSpec=','.join(frameSpec)
	for passNo in (1,2):
		cmd=['mencoder','mf://%s'%frameSpec,'-mf','fps=%d'%int(fps),'-ovc','lavc','-lavcopts','vbitrate=%d:vpass=%d:threads=%d:%s'%(int(kbps),passNo,O.numThreads,'turbo' if passNo==1 else ''),'-o',('/dev/null' if passNo==1 else out)]
		print('Pass %d:'%passNo,' '.join(cmd))
		ret=subprocess.call(cmd)
		if ret!=0: raise RuntimeError("Error when running mencoder.")

def replaceCollider(colliderEngine):
	"""Replaces collider (Collider) engine with the engine supplied. Raises error if no collider is in engines."""
	colliderIdx=-1
	for i,e in enumerate(O.engines):
		if O.isChildClassOf(e.__class__.__name__,"Collider"):
			colliderIdx=i
			break
	if colliderIdx<0: raise RuntimeError("No Collider found within O.engines.")
	O.engines=O.engines[:colliderIdx]+[colliderEngine]+O.engines[colliderIdx+1:]

def _procStatus(name):
	import os
	for l in open('/proc/%d/status'%os.getpid()):
		if l.split(':')[0]==name: return l
	raise("No such line in /proc/[pid]/status: "+name)
def vmData():
	"Return memory usage data from Linux's /proc/[pid]/status, line VmData."
	l=_procStatus('VmData'); ll=l.split(); assert(ll[2]=='kB')
	return int(ll[1])

def uniaxialTestFeatures(filename=None,areaSections=10,axis=-1,distFactor=2.2,**kw):
	"""Get some data about the current packing useful for uniaxial test:

#. Find the dimensions that is the longest (uniaxial loading axis)

#. Find the minimum cross-section area of the specimen by examining several (areaSections) sections perpendicular to axis, computing area of the convex hull for each one. This will work also for non-prismatic specimen.

#. Find the bodies that are on the negative/positive boundary, to which the straining condition should be applied.

:param filename: if given, disks will be loaded from this file (ASCII format); if not, current simulation will be used.
:param float areaSection: number of section that will be used to estimate cross-section
:param ∈{0,1,2} axis: if given, force strained axis, rather than computing it from predominant length
:return: dictionary with keys ``negIds``, ``posIds``, ``axis``, ``area``.

.. warning::
	The function :yref:`sudodem.utils.approxSectionArea` uses convex hull algorithm to find the area, but the implementation is reported to be *buggy* (bot works in some cases). Always check this number, or fix the convex hull algorithm (it is documented in the source, see :ysrc:`py/_utils.cpp`).

	"""
	if filename: ids=disksFromFile(filename,**kw)
	else: ids=[b.id for b in O.bodies]
	mm,mx=aabbExtrema()
	dim=aabbDim();
	if axis<0: axis=list(dim).index(max(dim)) # list(dim) for compat with python 2.5 which didn't have index defined for tuples yet (appeared in 2.6 first)
	assert(axis in (0,1,2))
	import numpy
	areas=[approxSectionArea(coord,axis) for coord in numpy.linspace(mm[axis],mx[axis],num=10)[1:-1]]
	negIds,posIds=negPosExtremeIds(axis=axis,distFactor=distFactor)
	return {'negIds':negIds,'posIds':posIds,'axis':axis,'area':min(areas)}

def voxelPorosityTriaxial(triax,resolution=200,offset=0):
	"""
	Calculate the porosity of a sample, given the TriaxialCompressionEngine.

	A function :yref:`sudodem.utils.voxelPorosity` is invoked, with the volume of a box enclosed by TriaxialCompressionEngine walls.
	The additional parameter offset allows using a smaller volume inside the box, where each side of the volume is at offset distance
	from the walls. By this way it is possible to find a more precise porosity of the sample, since at walls' contact the porosity is usually reduced.

	A recommended value of offset is bigger or equal to the average radius of disks inside.

	The value of resolution depends on size of disks used. It can be calibrated by invoking voxelPorosityTriaxial with offset=0 and
	comparing the result with TriaxialCompressionEngine.porosity. After calibration, the offset can be set to radius, or a bigger value, to get
	the result.

	:param triax: the TriaxialCompressionEngine handle
	:param resolution: voxel grid resolution
	:param offset: offset distance
	:return: the porosity of the sample inside given volume

	Example invocation::

		from sudodem import utils
		rAvg=0.03
		TriaxialTest(numberOfGrains=200,radiusMean=rAvg).load()
		O.dt=-1
		O.run(1000)
		O.engines[4].porosity
		0.44007807740143889
		utils.voxelPorosityTriaxial(O.engines[4],200,0)
		0.44055412500000002
		utils.voxelPorosityTriaxial(O.engines[4],200,rAvg)
		0.36798199999999998
	"""
	p_bottom	= O.bodies[triax.wall_bottom_id].state.se3[0]
	p_top		= O.bodies[triax.wall_top_id	 ].state.se3[0]
	p_left		= O.bodies[triax.wall_left_id	].state.se3[0]
	p_right		= O.bodies[triax.wall_right_id ].state.se3[0]
	p_front		= O.bodies[triax.wall_front_id ].state.se3[0]
	p_back		= O.bodies[triax.wall_back_id	].state.se3[0]
	th							= (triax.thickness)*0.5+offset
	x_0						 = p_left	[0] + th
	x_1						 = p_right [0] - th
	y_0						 = p_bottom[1] + th
	y_1						 = p_top	 [1] - th
	z_0						 = p_back	[2] + th
	z_1						 = p_front [2] - th
	a=Vector3(x_0,y_0,z_0)
	b=Vector3(x_1,y_1,z_1)
	return voxelPorosity(resolution,a,b)


def trackPerfomance(updateTime=5):
	"""
	Track perfomance of a simulation. (Experimental)
	Will create new thread to produce some plots.
	Useful for track perfomance of long run simulations (in bath mode for example).
	"""

	def __track_perfomance(updateTime):
		pid=os.getpid()
		threadsCpu={}
		lastTime,lastIter=-1,-1
		while 1:
			time.sleep(updateTime)
			if not O.running:
				lastTime,lastIter=-1,-1
				continue
			if lastTime==-1:
				lastTime=time.time();lastIter=O.iter
				plot.plots.update({'Iteration':('Perfomance',None,'Bodies','Interactions')})
				continue
			curTime=time.time();curIter=O.iter
			perf=(curIter-lastIter)/(curTime-lastTime)
			out=subprocess.Popen(['top','-bH','-n1', ''.join(['-p',str(pid)])],stdout=subprocess.PIPE).communicate()[0].splitlines()
			for s in out[7:-1]:
				w=s.split()
				threadsCpu[w[0]]=float(w[8])
			plot.addData(Iteration=curIter,Iter=curIter,Perfomance=perf,Bodies=len(O.bodies),Interactions=len(O.interactions),**threadsCpu)
			plot.plots.update({'Iter':threadsCpu.keys()})
			lastTime=time.time();lastIter=O.iter

	thread.start_new_thread(__track_perfomance,(updateTime))


def NormalRestitution2DampingRate(en):
	r"""Compute the normal damping rate as a function of the normal coefficient of restitution $e_n$. For $e_n\in\langle0,1\rangle$ damping rate equals

	.. math:: -\frac{\log e_n}{\sqrt{e_n^2+\pi^2}}

	"""
	if en == 0.0: return 0.999999999
	if en == 1.0: return 0.0
	from math import sqrt,log,pi
	ln_en = math.log(en)
	return (-ln_en/math.sqrt((math.pow(ln_en,2) + math.pi*math.pi)))

def xMirror(half):
	"""Mirror a sequence of 2d points around the x axis (changing sign on the y coord).
The sequence should start up and then it will wrap from y downwards (or vice versa).
If the last point's x coord is zero, it will not be duplicated."""
	return list(half)+[(x,-y) for x,y in reversed(half[:-1] if half[-1][1]==0 else half)]

#############################
##### deprecated functions


def _deprecatedUtilsFunction(old,new):
	"Wrapper for deprecated functions, example below."
	import warnings
	warnings.warn('Function utils.%s is deprecated, use %s instead.'%(old,new),stacklevel=2,category=DeprecationWarning)

# example of _deprecatedUtilsFunction usage:
#
# def import_mesh_geometry(*args,**kw):
#		"|ydeprecated|"
#		_deprecatedUtilsFunction('import_mesh_geometry','sudodem.import.gmsh')
#		import sudodem.ymport
#		return sudodem.ymport.stl(*args,**kw)


class TableParamReader():
	"""Class for reading simulation parameters from text file.

Each parameter is represented by one column, each parameter set by one line. Colums are separated by blanks (no quoting).

First non-empty line contains column titles (without quotes).
You may use special column named 'description' to describe this parameter set;
if such colum is absent, description will be built by concatenating column names and corresponding values (``param1=34,param2=12.22,param4=foo``)

* from columns ending in ``!`` (the ``!`` is not included in the column name)
* from all columns, if no columns end in ``!``.

Empty lines within the file are ignored (although counted); ``#`` starts comment till the end of line. Number of blank-separated columns must be the same for all non-empty lines.

A special value ``=`` can be used instead of parameter value; value from the previous non-empty line will be used instead (works recursively).

This class is used by :yref:`sudodem.utils.readParamsFromTable`.
	"""
	def __init__(self,file):
		"Setup the reader class, read data into memory."
		import re
		# read file in memory, remove newlines and comments; the [''] makes lines 1-indexed
		ll=[re.sub('\s*#.*','',l[:-1]) for l in ['']+open(file,'r').readlines()]
		# usable lines are those that contain something else than just spaces
		usableLines=[i for i in range(len(ll)) if not re.match(r'^\s*(#.*)?$',ll[i])]
		headings=ll[usableLines[0]].split()
		# use all values of which heading has ! after its name to build up the description string
		# if there are none, use all columns
		if not 'description' in headings:
			bangHeads=[h[:-1] for h in headings if h[-1]=='!'] or headings
			headings=[(h[:-1] if h[-1]=='!' else h) for h in headings]
		usableLines=usableLines[1:] # and remove headinds from usableLines
		values={}
		for l in usableLines:
			val={}
			for i in range(len(headings)):
				val[headings[i]]=ll[l].split()[i]
			values[l]=val
		lines=values.keys(); lines.sort()
		# replace '=' by previous value of the parameter
		for i,l in enumerate(lines):
			for j in values[l].keys():
				if values[l][j]=='=':
					try:
						values[l][j]=values[lines[i-1]][j]
					except(IndexError,KeyError):
						raise RuntimeError("The = specifier on line %d refers to nonexistent value on previous line?"%l)
		#import pprint; pprint.pprint(headings); pprint.pprint(values)
		# add descriptions, but if they repeat, append line number as well
		if not 'description' in headings:
			descs=set()
			for l in lines:
				dd=','.join(head.replace('!','')+'='+('%g'%values[head] if isinstance(values[l][head],float) else str(values[l][head])) for head in bangHeads).replace("'",'').replace('"','')
				if dd in descs: dd+='__line=%d__'%l
				values[l]['description']=dd
				descs.add(dd)
		self.values=values

	def paramDict(self):
		"""Return dictionary containing data from file given to constructor. Keys are line numbers (which might be non-contiguous and refer to real line numbers that one can see in text editors), values are dictionaries mapping parameter names to their values given in the file. The special value '=' has already been interpreted, ``!`` (bangs) (if any) were already removed from column titles, ``description`` column has already been added (if absent)."""
		return self.values

if __name__=="__main__":
	tryTable="""head1 important2! !OMP_NUM_THREADS! abcd
	1 1.1 1.2 1.3
	'a' 'b' 'c' 'd'	### comment

	# empty line
	1 = = g
"""
	file='/tmp/try-tbl.txt'
	f=open(file,'w')
	f.write(tryTable)
	f.close()
	from pprint import *
	pprint(TableParamReader(file).paramDict())

def runningInBatch():
	'Tell whether we are running inside the batch or separately.'
	import os
	return 'SUDODEM_BATCH' in os.environ

def waitIfBatch():
	'Block the simulation if running inside a batch. Typically used at the end of script so that it does not finish prematurely in batch mode (the execution would be ended in such a case).'
	if runningInBatch(): O.wait()

def readParamsFromTable(tableFileLine=None,noTableOk=True,unknownOk=False,**kw):
	"""
	Read parameters from a file and assign them to __builtin__ variables.

	The format of the file is as follows (commens starting with # and empty lines allowed)::

		# commented lines allowed anywhere
		name1 name2 … # first non-blank line are column headings
					# empty line is OK, with or without comment
		val1	val2	… # 1st parameter set
		val2	val2	… # 2nd
		…

	Assigned tags (the ``description`` column is synthesized if absent,see :yref:`sudodem.utils.TableParamReader`);

		O.tags['description']=…																			# assigns the description column; might be synthesized
		O.tags['params']="name1=val1,name2=val2,…"									 # all explicitly assigned parameters
		O.tags['defaultParams']="unassignedName1=defaultValue1,…"		# parameters that were left at their defaults
		O.tags['d.id']=O.tags['id']+'.'+O.tags['description']
		O.tags['id.d']=O.tags['description']+'.'+O.tags['id']

	All parameters (default as well as settable) are saved using :yref:`sudodem.utils.saveVars`\ ``('table')``.

	:param tableFile: text file (with one value per blank-separated columns)
	:param int tableLine: number of line where to get the values from
	:param bool noTableOk: if False, raise exception if the file cannot be open; use default values otherwise
	:param bool unknownOk: do not raise exception if unknown column name is found in the file, and assign it as well
	:return: number of assigned parameters
	"""
	tagsParams=[]
	# dictParams is what eventually ends up in sudodem.params.table (default+specified values)
	dictDefaults,dictParams,dictAssign={},{},{}
	import os, __builtin__,re,math
	if not tableFileLine and ('SUDODEM_BATCH' not in os.environ or os.environ['SUDODEM_BATCH']==''):
		if not noTableOk: raise EnvironmentError("SUDODEM_BATCH is not defined in the environment")
		O.tags['line']='l!'
	else:
		if not tableFileLine: tableFileLine=os.environ['SUDODEM_BATCH']
		env=tableFileLine.split(':')
		tableFile,tableLine=env[0],int(env[1])
		allTab=TableParamReader(tableFile).paramDict()
		if not allTab.has_key(tableLine): raise RuntimeError("Table %s doesn't contain valid line number %d"%(tableFile,tableLine))
		vv=allTab[tableLine]
		O.tags['line']='l%d'%tableLine
		O.tags['description']=vv['description']
		O.tags['id.d']=O.tags['id']+'.'+O.tags['description']; O.tags['d.id']=O.tags['description']+'.'+O.tags['id']
		# assign values specified in the table to python vars
		# !something cols are skipped, those are env vars we don't treat at all (they are contained in description, though)
		for col in vv.keys():
			if col=='description' or col[0]=='!': continue
			if col not in kw.keys() and (not unknownOk): raise NameError("Parameter `%s' has no default value assigned"%col)
			if vv[col]=='*': vv[col]=kw[col] # use default value for * in the table
			elif col in kw.keys(): kw.pop(col) # remove the var from kw, so that it contains only those that were default at the end of this loop
			#print 'ASSIGN',col,vv[col]
			tagsParams+=['%s=%s'%(col,vv[col])];
			dictParams[col]=eval(vv[col],math.__dict__)
	# assign remaining (default) keys to python vars
	defaults=[]
	for k in kw.keys():
		dictDefaults[k]=kw[k]
		defaults+=["%s=%s"%(k,kw[k])];
	O.tags['defaultParams']=",".join(defaults)
	O.tags['params']=",".join(tagsParams)
	dictParams.update(dictDefaults)
	saveVars('table',loadNow=True,**dictParams)
	return len(tagsParams)

def psd(bins=5, mass=True, mask=-1):
	"""Calculates particle size distribution.

	:param int bins: number of bins
	:param bool mass: if true, the mass-PSD will be calculated
	:param int mask: :yref:`Body.mask` for the body
	:return:
		* binsSizes: list of bin's sizes
		* binsProc: how much material (in percents) are in the bin, cumulative
		* binsSumCum: how much material (in units) are in the bin, cumulative

		binsSizes, binsProc, binsSumCum
	"""
	maxD = 0.0
	minD = 0.0

	for b in O.bodies:
		if (isinstance(b.shape,Disk) and ((mask<0) or ((b.mask&mask)!=0))):
			if ((2*b.shape.radius)	> maxD) : maxD = 2*b.shape.radius
			if (((2*b.shape.radius)	< minD) or (minD==0.0)): minD = 2*b.shape.radius

	if (minD==maxD): return false       #All particles are having the same size

	binsSizes = numpy.linspace(minD, maxD, bins+1)

	deltaBinD = (maxD-minD)/bins
	binsMass = numpy.zeros(bins)
	binsNumbers = numpy.zeros(bins)

	for b in O.bodies:
		if (isinstance(b.shape,Disk) and ((mask<0) or ((b.mask&mask)!=0))):
			d=2*b.shape.radius

			basketId = int(math.floor( (d-minD) / deltaBinD ) )
			if (d == maxD): basketId = bins-1												 #If the diameter equals the maximal diameter, put the particle into the last bin
			binsMass[basketId] = binsMass[basketId] + b.state.mass		#Put masses into the bin
			binsNumbers[basketId] = binsNumbers[basketId] + 1					#Put numbers into the bin

	binsProc = numpy.zeros(bins+1)
	binsSumCum = numpy.zeros(bins+1)

	i=1
	for size in binsSizes[:-1]:
		if (mass): binsSumCum[i]+=binsSumCum[i-1]+binsMass[i-1]		#Calculate mass
		else: binsSumCum[i]+=binsSumCum[i-1]+binsNumbers[i-1]			#Calculate number of particles
		i+=1

	if (binsSumCum[len(binsSumCum)-1] > 0):
		i=0
		for l in binsSumCum:
			binsProc[i] = binsSumCum[i]/binsSumCum[len(binsSumCum)-1]
			i+=1
	return binsSizes, binsProc, binsSumCum

class clumpTemplate:
	"""Create a clump template by a list of relative radii and a list of relative positions. Both lists must have the same length.

	:param [float,float,...] relRadii: list of relative radii (minimum length = 2)
	:param [Vector3,Vector3,...] relPositions: list of relative positions (minimum length = 2)

	"""
	def __init__(self,relRadii=[],relPositions=[[],[]]):
		if (len(relRadii) != len(relPositions)):
			raise ValueError("Given lists must have same length! Given lists does not match clump template structure.");
		if (len(relRadii) < 2):
			raise ValueError("One or more of given lists for relative radii have length < 2! Given lists does not match clump template structure.");
		for ii in range(0,len(relPositions)):
			if len(relPositions[ii]) != 3:
				raise ValueError("One or more of given lists for relative positions do not have length of 3! Given lists does not match clump template structure.");
			for jj in range(ii,len(relPositions)):
				if ii != jj:
					if (relPositions[ii] == relPositions[jj]):
						raise ValueError("Two or more of given lists for relative positions are equal! Given lists does not match clump template structure.");
		self.numCM = len(relRadii)
		self.relRadii = relRadii
		self.relPositions = relPositions

class UnstructuredGrid:
	"""EXPERIMENTAL.
	Class representing triangulated FEM-like unstructured grid. It is used for transfereing data from ad to SUDODEM and external FEM program. The main purpose of this class is to store information about individual grid vertices/nodes coords (since facets stores only coordinates of vertices in local coords) and to avaluate and/or apply nodal forces from contact forces (from actual contact force and contact point the force is distributed to nodes using linear approximation).

	TODO rewrite to C++
	TODO better docs

	:param dict vertices: dict of {internal vertex label:vertex}, e.g. {5:(0,0,0),22:(0,1,0),23:(1,0,0)}
	:param dict connectivityTable: dict of {internal element label:[indices of vertices]}, e.g. {88:[5,22,23]}
	"""
	def __init__(self,vertices={},connectivityTable={},**kw):
		self.vertices = vertices
		self.connectivityTable = connectivityTable
		self.forces = dict( (i,Vector3(0,0,0)) for i in vertices)
		self.build(**kw)
	def setup(self,vertices,connectivityTable,toSimulation=False,bodies=None,**kw):
		"""Sets new information to receiver

		:param dict vertices: see constructor for explanation
		:param dict connectivityTable: see constructor for explanation
		:param bool toSimulation: if new information should be inserted to SudoDEM simulation (create new bodies or not)
		:param [[int]]|None bodies: list of list of bodies indices to be appended as clumps (thus no contact detection is done within one body)
		"""
		self.vertices = dict( (i,v) for i,v in vertices.iteritems())
		self.connectivityTable = dict( (i,e) for i,e in connectivityTable.iteritems())
		self.build(**kw)
		if toSimulation:
			self.toSimulation(bodies)
	def build(self,**kw):
		self.elements = {}
		for i,c in self.connectivityTable.iteritems():
			if len(c) == 3:
				b = facet([self.vertices[j] for j in c],**kw)
			elif len(c) == 4:
				b = tetra([self.vertices[j] for j in c],**kw)
				#b = polyhedron([self.vertices[j] for j in c],**kw)
			else:
				raise(RuntimeError, "Unsupported cell shape (should be triangle or tetrahedron)")
			self.elements[i] = b
	def resetForces(self):
		for i in self.vertices:
			self.forces[i] = Vector3(0,0,0)
	def getForcesOfNodes(self):
		"""Computes forces for each vertex/node. The nodal force is computed from contact force and contact point using linear approximation
		"""
		self.resetForces()
		for i,e in self.elements.iteritems():
			ie = self.connectivityTable[i]
			for i in e.intrs():
				fn = i.phys.normalForce
				fs = i.phys.shearForce if hasattr(i.phys,"shearForce") else Vector3.Zero
				f = (fn+fs) * (-1 if e.id == i.id1 else +1 if e.id == i.id2 else 'ERROR')
				cp = i.geom.contactPoint
				if isinstance(e.shape,Facet):
					v0,v1,v2 = [Vector3(self.vertices[j]) for j in ie]
					w0 = ((cp-v1).cross(v2-v1)).norm()
					w1 = ((cp-v2).cross(v0-v2)).norm()
					w2 = ((cp-v0).cross(v1-v0)).norm()
					ww = w0+w1+w2
					self.forces[ie[0]] += f*w0/ww
					self.forces[ie[1]] += f*w1/ww
					self.forces[ie[2]] += f*w2/ww
				elif isinstance(e.shape,Tetra):
					v0,v1,v2,v3 = [Vector3(self.vertices[j]) for j in ie]
					w0 = TetrahedronVolume((cp,v1,v2,v3))
					w1 = TetrahedronVolume((cp,v2,v3,v0))
					w2 = TetrahedronVolume((cp,v3,v0,v1))
					w3 = TetrahedronVolume((cp,v0,v1,v2))
					ww = w0+w1+w2+w3
					self.forces[ie[0]] += f*w0/ww
					self.forces[ie[1]] += f*w1/ww
					self.forces[ie[2]] += f*w2/ww
					self.forces[ie[3]] += f*w3/ww
				else:
					raise(RuntimeError, "TODO")
		return self.forces
	def setPositionsOfNodes(self,newPoss):
		"""Sets new position of nodes and also updates all elements in the simulation

		:param [Vector3] newPoss: list of new positions
		"""
		for i in self.vertices:
			self.vertices[i] = newPoss[i]
		self.updateElements()
	def updateElements(self):
		"""Updates positions of all elements in the simulation
		"""
		for i,c in self.connectivityTable.iteritems():
			e = self.elements[i]
			e.state.pos = Vector3(0,0,0)
			e.state.ori = Quaternion((1,0,0),0)
			if isinstance(e.shape,Facet):
				e.shape.vertices = [self.vertices[j] for j in c]
			elif isinstance(e.shape,Tetra):
				e.shape.v = [self.vertices[j] for j in c]
			else:
				raise(RuntimeError, "TODO")
	def toSimulation(self,bodies=None):
		"""Insert all elements to SudoDEM simulation
		"""
		if bodies:
			e = self.elements.values()
			for b in bodies:
				O.bodies.appendClumped([e[i] for i in b])
		else:
			O.bodies.append(self.elements.values())
