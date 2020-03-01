# encoding: utf-8
# 2009 © Václav Šmilauer <eudoxos@arcig.cz>

"""
Creating packings and filling volumes defined by boundary representation or constructive solid geometry.

For examples, see

* :ysrc:`scripts/test/gts-operators.py`
* :ysrc:`scripts/test/gts-random-pack-obb.py`
* :ysrc:`scripts/test/gts-random-pack.py`
* :ysrc:`scripts/test/pack-cloud.py`
* :ysrc:`scripts/test/pack-predicates.py`
* :ysrc:`examples/packs/packs.py`
* :ysrc:`examples/gts-horse/gts-horse.py`
* :ysrc:`examples/WireMatPM/wirepackings.py`
"""

import itertools,warnings
from numpy import arange
from math import sqrt
from sudodem import utils

from sudodem.wrapper import *

try:
	from minieigen import *
except ImportError:
	from miniEigen import *

## compatibility hack for python 2.5 (21/8/2009)
## can be safely removed at some point
if 'product' not in dir(itertools):
	def product(*args, **kwds):
		"http://docs.python.org/library/itertools.html#itertools.product"
		pools = map(tuple, args) * kwds.get('repeat', 1); result = [[]]
		for pool in pools: result = [x+[y] for x in result for y in pool]
		for prod in result: yield tuple(prod)
	itertools.product=product

# for now skip the import, but try in inGtsSurface constructor again, to give error if we really use it
try:
	import gts
except ImportError: pass

# make c++ predicates available in this module
noPredicate = False
try:
	from _packPredicates import * ## imported in randomDensePack as well
except ImportError: pass; noPredicate = True

# import DiskPack
#from _packDisks import *
#from _packObb import *

##
# extend _packDisk.DiskPack c++ class by this method
##
if not (noPredicate):
  def DiskPack_toSimulation(self,rot=Matrix3.Identity,**kw):
    """Append disks directly to the simulation. In addition calling :yref:`O.bodies.append<BodyContainer.append>`,
  this method also appropriately sets periodic cell information of the simulation.

    >>> from sudodem import pack; from math import *
    >>> sp=pack.DiskPack()

  Create random periodic packing with 20 disks:

    >>> sp.makeCloud((0,0,0),(5,5,5),rMean=.5,rRelFuzz=.5,periodic=True,num=20)
    20

  Virgin simulation is aperiodic:

    >>> O.reset()
    >>> O.periodic
    False

  Add generated packing to the simulation, rotated by 45° along +z

    >>> sp.toSimulation(rot=Quaternion((0,0,1),pi/4),color=(0,0,1))
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

  Periodic properties are transferred to the simulation correctly, including rotation (this could be avoided by explicitly passing "hSize=O.cell.hSize" as an argument):

    >>> O.periodic
    True
    >>> O.cell.refSize
    Vector3(5,5,5)
    """
  #The following 2 lines do not work, because of accuaracy
  #>>> O.cell.hSize
  #Matrix3(3.53553,-3.53553,0, 3.53553,3.53553,0, 0,0,5)
    """
  The current state (even if rotated) is taken as mechanically undeformed, i.e. with identity transformation:

    >>> O.cell.trsf
    Matrix3(1,0,0, 0,1,0, 0,0,1)

  :param Quaternion/Matrix3 rot: rotation of the packing, which will be applied on disks and will be used to set :yref:`Cell.trsf` as well.
  :param \*\*kw: passed to :yref:`sudodem.utils.disk`
  :return: list of body ids added (like :yref:`O.bodies.append<BodyContainer.append>`)
  """
    if isinstance(rot,Quaternion): rot=rot.toRotationMatrix()
    assert(isinstance(rot,Matrix3))
    if self.isPeriodic: O.periodic=True
    if self.cellSize!=Vector3.Zero and self.isPeriodic:
      O.cell.hSize=rot*Matrix3(self.cellSize[0],0,0, 0,self.cellSize[1],0, 0,0,self.cellSize[2])
      O.cell.trsf=Matrix3.Identity
    if not self.hasClumps():
      return O.bodies.append([utils.disk(rot*c,r,**kw) for c,r in self])
    else:
      standalone,clumps=self.getClumps()
      ids=O.bodies.append([utils.disk(rot*c,r,**kw) for c,r in self]) # append all disks first
      clumpIds=[]
      userColor='color' in kw
      for clump in clumps:
        clumpIds.append(O.bodies.clump([ids[i] for i in clump])) # clump disks with given ids together, creating the clump object as well
        # make all disks within one clump a single color, unless color was specified by the user
        if not userColor:
          for i in clump[1:]: O.bodies[ids[i]].shape.color=O.bodies[ids[clump[0]]].shape.color
      return ids+clumpIds

  DiskPack.toSimulation=DiskPack_toSimulation


  class inGtsSurface_py(Predicate):
    """This class was re-implemented in c++, but should stay here to serve as reference for implementing
    Predicates in pure python code. C++ allows us to play dirty tricks in GTS which are not accessible
    through pygts itself; the performance penalty of pygts comes from fact that if constructs and destructs
    bb tree for the surface at every invocation of gts.Point().is_inside(). That is cached in the c++ code,
    provided that the surface is not manipulated with during lifetime of the object (user's responsibility).

    ---

    Predicate for GTS surfaces. Constructed using an already existing surfaces, which must be closed.

      import gts
      surf=gts.read(open('horse.gts'))
      inGtsSurface(surf)

    .. note::
      Padding is optionally supported by testing 6 points along the axes in the pad distance. This
      must be enabled in the ctor by saying doSlowPad=True. If it is not enabled and pad is not zero,
      warning is issued.
    """
    def __init__(self,surf,noPad=False):
      # call base class ctor; necessary for virtual methods to work as expected.
      # see comments in _packPredicates.cpp for struct PredicateWrap.
      super(inGtsSurface,self).__init__()
      if not surf.is_closed(): raise RuntimeError("Surface for inGtsSurface predicate must be closed.")
      self.surf=surf
      self.noPad=noPad
      inf=float('inf')
      mn,mx=[inf,inf,inf],[-inf,-inf,-inf]
      for v in surf.vertices():
        c=v.coords()
        mn,mx=[min(mn[i],c[i]) for i in 0,1,2],[max(mx[i],c[i]) for i in 0,1,2]
      self.mn,self.mx=tuple(mn),tuple(mx)
      import gts
    def aabb(self): return self.mn,self.mx
    def __call__(self,_pt,pad=0.):
      p=gts.Point(*_pt)
      if self.noPad:
        if pad!=0: warnings.warn("Padding disabled in ctor, using 0 instead.")
        return p.is_inside(self.surf)
      pp=[gts.Point(_pt[0]-pad,_pt[1],_pt[2]),gts.Point(_pt[0]+pad,_pt[1],_pt[2]),gts.Point(_pt[0],_pt[1]-pad,_pt[2]),gts.Point(_pt[0],_pt[1]+pad,_pt[2]),gts.Point(_pt[0],_pt[1],_pt[2]-pad),gts.Point(_pt[0],_pt[1],_pt[2]+pad)]
      return p.is_inside(self.surf) and pp[0].is_inside(self.surf) and pp[1].is_inside(self.surf) and pp[2].is_inside(self.surf) and pp[3].is_inside(self.surf) and pp[4].is_inside(self.surf) and pp[5].is_inside(self.surf)

  class inSpace(Predicate):
    """Predicate returning True for any points, with infinite bounding box."""
    def __init__(self, _center=Vector3().Zero): self._center=_center
    def aabb(self):
      inf=float('inf'); return Vector3(-inf,-inf,-inf),Vector3(inf,inf,inf)
    def center(self): return self._center
    def dim(self):
      inf=float('inf'); return Vector3(inf,inf,inf)
    def __call__(self,pt,pad): return True

#####
## surface construction and manipulation
#####

def gtsSurface2Facets(surf,**kw):
	"""Construct facets from given GTS surface. \*\*kw is passed to utils.facet."""
	import gts
	return [utils.facet([v.coords() for v in face.vertices()],**kw) for face in surf.faces()]

def sweptPolylines2gtsSurface(pts,threshold=0,capStart=False,capEnd=False):
	"""Create swept suface (as GTS triangulation) given same-length sequences of points (as 3-tuples).

If threshold is given (>0), then

* degenerate faces (with edges shorter than threshold) will not be created
* gts.Surface().cleanup(threshold) will be called before returning, which merges vertices mutually closer than threshold. In case your pts are closed (last point concident with the first one) this will the surface strip of triangles. If you additionally have capStart==True and capEnd==True, the surface will be closed.

.. note:: capStart and capEnd make the most naive polygon triangulation (diagonals) and will perhaps fail for non-convex sections.

.. warning:: the algorithm connects points sequentially; if two polylines are mutually rotated or have inverse sense, the algorithm will not detect it and connect them regardless in their given order.
	"""
	import gts # will raise an exception in gts-less builds
	if not len(set([len(pts1) for pts1 in pts]))==1: raise RuntimeError("Polylines must be all of the same length!")
	vtxs=[[gts.Vertex(x,y,z) for x,y,z in pts1] for pts1 in pts]
	sectEdges=[[gts.Edge(vtx[i],vtx[i+1]) for i in xrange(0,len(vtx)-1)] for vtx in vtxs]
	interSectEdges=[[] for i in range(0,len(vtxs)-1)]
	for i in range(0,len(vtxs)-1):
		for j in range(0,len(vtxs[i])):
			interSectEdges[i].append(gts.Edge(vtxs[i][j],vtxs[i+1][j]))
			if j<len(vtxs[i])-1: interSectEdges[i].append(gts.Edge(vtxs[i][j],vtxs[i+1][j+1]))
	if threshold>0: # replace edges of zero length with None; their faces will be skipped
		def fixEdges(edges):
			for i,e in enumerate(edges):
				if (Vector3(e.v1.x,e.v1.y,e.v1.z)-Vector3(e.v2.x,e.v2.y,e.v2.z)).norm()<threshold: edges[i]=None
		for ee in sectEdges: fixEdges(ee)
		for ee in interSectEdges: fixEdges(ee)
	surf=gts.Surface()
	for i in range(0,len(vtxs)-1):
		for j in range(0,len(vtxs[i])-1):
			ee1=interSectEdges[i][2*j+1],sectEdges[i+1][j],interSectEdges[i][2*j]
			ee2=sectEdges[i][j],interSectEdges[i][2*j+2],interSectEdges[i][2*j+1]
			if None not in ee1: surf.add(gts.Face(interSectEdges[i][2*j+1],sectEdges[i+1][j],interSectEdges[i][2*j]))
			if None not in ee2: surf.add(gts.Face(sectEdges[i][j],interSectEdges[i][2*j+2],interSectEdges[i][2*j+1]))
	def doCap(vtx,edg,start):
		ret=[]
		eFan=[edg[0]]+[gts.Edge(vtx[i],vtx[0]) for i in range(2,len(vtx))]
		for i in range(1,len(edg)):
			ret+=[gts.Face(eFan[i-1],eFan[i],edg[i]) if start else gts.Face(eFan[i-1],edg[i],eFan[i])]
		return ret
	caps=[]
	if capStart: caps+=doCap(vtxs[0],sectEdges[0],start=True)
	if capEnd: caps+=doCap(vtxs[-1],sectEdges[-1],start=False)
	for cap in caps: surf.add(cap)
	if threshold>0: surf.cleanup(threshold)
	return surf

def gtsSurfaceBestFitOBB(surf):
	"""Return (Vector3 center, Vector3 halfSize, Quaternion orientation) describing
	best-fit oriented bounding box (OBB) for the given surface. See cloudBestFitOBB
	for details."""
	import gts
	pts=[Vector3(v.x,v.y,v.z) for v in surf.vertices()]
	return cloudBestFitOBB(tuple(pts))

def revolutionSurfaceMeridians(sects,angles,origin=Vector3().Zero,orientation=Quaternion().Identity):
	"""Revolution surface given sequences of 2d points and sequence of corresponding angles,
	returning sequences of 3d points representing meridian sections of the revolution surface.
	The 2d sections are turned around z-axis, but they can be transformed
	using the origin and orientation arguments to give arbitrary orientation."""
	import math
	def toGlobal(x,y,z):
		return tuple(origin+orientation*(Vector3(x,y,z)))
	return [[toGlobal(x2d*math.cos(angles[i]),x2d*math.sin(angles[i]),y2d) for x2d,y2d in sects[i]] for i in range(0,len(sects))]

########
## packing generators
########


def regularOrtho(predicate,radius,gap,**kw):
	"""Return set of disks in regular orthogonal grid, clipped inside solid given by predicate.
	Created disks will have given radius and will be separated by gap space."""
	ret=[]
	mn,mx=predicate.aabb()
	if(max([mx[i]-mn[i] for i in 0,1,2])==float('inf')): raise ValueError("Aabb of the predicate must not be infinite (didn't you use union | instead of intersection & for unbounded predicate such as notInNotch?");
	xx,yy,zz=[arange(mn[i]+radius,mx[i]-radius,2*radius+gap) for i in 0,1,2]
	for xyz in itertools.product(xx,yy,zz):
		if predicate(xyz,radius): ret+=[utils.disk(xyz,radius=radius,**kw)]
	if (len(ret)==0):
		warnings.warn('No disks are produced by regularOrtho-function',category=RuntimeWarning)
	return ret

def regularHexa(predicate,radius,gap,**kw):
	"""Return set of disks in regular hexagonal grid, clipped inside solid given by predicate.
	Created disks will have given radius and will be separated by gap space."""
	ret=[]
	a=2*radius+gap
	# thanks to Nasibeh Moradi for finding bug here:
	# http://www.mail-archive.com/sudodem-users@lists.launchpad.net/msg01424.html
	hy,hz=a*sqrt(3)/2.,a*sqrt(6)/3.
	mn,mx=predicate.aabb()
	dim=[mx[i]-mn[i] for i in 0,1,2]
	if(max(dim)==float('inf')): raise ValueError("Aabb of the predicate must not be infinite (didn't you use union | instead of intersection & for unbounded predicate such as notInNotch?");
	ii,jj,kk=[range(0,int(dim[0]/a)+1),range(0,int(dim[1]/hy)+1),range(0,int(dim[2]/hz)+1)]
	for i,j,k in itertools.product(ii,jj,kk):
		x,y,z=mn[0]+radius+i*a,mn[1]+radius+j*hy,mn[2]+radius+k*hz
		if j%2==0: x+= a/2. if k%2==0 else -a/2.
		if k%2!=0: x+=a/2.; y+=hy/2.
		if predicate((x,y,z),radius): ret+=[utils.disk((x,y,z),radius=radius,**kw)]
	if (len(ret)==0):
		warnings.warn('No disks are produced by regularHexa-function',category=RuntimeWarning)
	return ret

def filterDiskPack(predicate,diskPack,returnDiskPack=None,**kw):
	"""Using given DiskPack instance, return disks that satisfy predicate.
	It returns either a :yref:`sudodem._packDisks.DiskPack` (if returnDiskPack) or a list.
	The packing will be recentered to match the predicate and warning is given if the predicate
	is larger than the packing."""
	if returnDiskPack==None:
		warnings.warn('The default behavior will change; specify returnDiskPack=True for the new behavior, and False to get rid of this warning (your code will break in the future, however). The returned DiskPack object can be added to the simulation using DiskPack.toSimulation()',category=FutureWarning)
		returnDiskPack=False
	mn,mx=predicate.aabb()
	dimP,centP=predicate.dim(),predicate.center()
	dimS,centS=diskPack.dim(),diskPack.center()
	if dimP[0]>dimS[0] or dimP[1]>dimS[1] or dimP[2]>dimS[2]: warnings.warn("Packing's dimension (%s) doesn't fully contain dimension of the predicate (%s)."%(dimS,dimP))
	diskPack.translate(centP-centS)
	if returnDiskPack:
		ret=DiskPack()
		for c,r in diskPack:
			if predicate(c,r): ret.add(c,r)
		return ret
	else:
		# return particles to be added to O.bodies
		ret=[]
		for s in diskPack:
			if predicate(s[0],s[1]): ret+=[utils.disk(s[0],radius=s[1],**kw)]
		return ret

def _memoizePacking(memoizeDb,sp,radius,rRelFuzz,wantPeri,fullDim,noPrint=False):
	if not memoizeDb: return
	import cPickle,sqlite3,time,os
	if os.path.exists(memoizeDb):
		conn=sqlite3.connect(memoizeDb)
	else:
		conn=sqlite3.connect(memoizeDb)
		c=conn.cursor()
		c.execute('create table packings (radius real, rRelFuzz real, dimx real, dimy real, dimz real, N integer, timestamp real, periodic integer, pack blob)')
	c=conn.cursor()
	packBlob=buffer(cPickle.dumps(sp.toList(),cPickle.HIGHEST_PROTOCOL))
	packDim=sp.cellSize if wantPeri else fullDim
	c.execute('insert into packings values (?,?,?,?,?,?,?,?,?)',(radius,rRelFuzz,packDim[0],packDim[1],packDim[2],len(sp),time.time(),wantPeri,packBlob,))
	c.close()
	conn.commit()
	if not noPrint: print "Packing saved to the database",memoizeDb

def _getMemoizedPacking(memoizeDb,radius,rRelFuzz,x1,y1,z1,fullDim,wantPeri,fillPeriodic,disksInCell,memoDbg=False,noPrint=False):
	"""Return suitable DiskPack read from *memoizeDb* if found, None otherwise.

		:param fillPeriodic: whether to fill fullDim by repeating periodic packing
		:param wantPeri: only consider periodic packings
	"""
	import os,os.path,sqlite3,time,cPickle,sys
	if memoDbg and not noPrint:
		def memoDbgMsg(s): print s
	else:
		def memoDbgMsg(s): pass
	if not memoizeDb or not os.path.exists(memoizeDb):
		if memoizeDb: memoDbgMsg("Database %s does not exist."%memoizeDb)
		return None
	# find suitable packing and return it directly
	conn=sqlite3.connect(memoizeDb); c=conn.cursor();
	try:
		c.execute('select radius,rRelFuzz,dimx,dimy,dimz,N,timestamp,periodic from packings order by N')
	except sqlite3.OperationalError:
		raise RuntimeError("ERROR: database `"+memoizeDb+"' not compatible with randomDensePack (corrupt, deprecated format or not a db created by randomDensePack)")
	for row in c:
		R,rDev,X,Y,Z,NN,timestamp,isPeri=row[0:8]; scale=radius/R
		rDev*=scale; X*=scale; Y*=scale; Z*=scale
		memoDbgMsg("Considering packing (radius=%g±%g,N=%g,dim=%g×%g×%g,%s,scale=%g), created %s"%(R,.5*rDev,NN,X,Y,Z,"periodic" if isPeri else "non-periodic",scale,time.asctime(time.gmtime(timestamp))))
		if not isPeri and wantPeri: memoDbgMsg("REJECT: is not periodic, which is requested."); continue
		if wantPeri and (X/x1>0.9 or X/x1<0.6):  memoDbgMsg("REJECT: initSize differs too much from scaled packing size."); continue
		if (rRelFuzz==0 and rDev!=0) or (rRelFuzz!=0 and rDev==0) or (rRelFuzz!=0 and abs((rDev-rRelFuzz)/rRelFuzz)>1e-2): memoDbgMsg("REJECT: radius fuzz differs too much (%g, %g desired)"%(rDev,rRelFuzz)); continue # radius fuzz differs too much
		if isPeri and wantPeri:
			if disksInCell>NN and disksInCell>0: memoDbgMsg("REJECT: Number of disks in the packing too small"); continue
			if abs((y1/x1)/(Y/X)-1)>0.3 or abs((z1/x1)/(Z/X)-1)>0.3: memoDbgMsg("REJECT: proportions (y/x=%g, z/x=%g) differ too much from what is desired (%g, %g)."%(Y/X,Z/X,y1/x1,z1/x1)); continue
		else:
			if (X<fullDim[0] or Y<fullDim[1] or Z<fullDim[2]): memoDbgMsg("REJECT: not large enough"); continue # not large enough
		memoDbgMsg("ACCEPTED");
		if not noPrint: print "Found suitable packing in %s (radius=%g±%g,N=%g,dim=%g×%g×%g,%s,scale=%g), created %s"%(memoizeDb,R,rDev,NN,X,Y,Z,"periodic" if isPeri else "non-periodic",scale,time.asctime(time.gmtime(timestamp)))
		c.execute('select pack from packings where timestamp=?',(timestamp,))
		sp=DiskPack(cPickle.loads(str(c.fetchone()[0])))
		sp.scale(scale);
		if isPeri and wantPeri:
			sp.isPeriodic = True
			sp.cellSize=(X,Y,Z);
			if fillPeriodic: sp.cellFill(Vector3(fullDim[0],fullDim[1],fullDim[2]));
		#sp.cellSize=(0,0,0) # resetting cellSize avoids warning when rotating
		return sp
		#if orientation: sp.rotate(*orientation.toAxisAngle())
		#return filterDiskPack(predicate,sp,material=material)
	#print "No suitable packing in database found, running",'PERIODIC compression' if wantPeri else 'triaxial'
	#sys.stdout.flush()



def randomDensePack(predicate,radius,material=-1,dim=None,cropLayers=0,rRelFuzz=0.,disksInCell=0,memoizeDb=None,useOBB=False,memoDbg=False,color=None,returnDiskPack=None):
	"""Generator of random dense packing with given geometry properties, using TriaxialTest (aperiodic)
	or PeriIsoCompressor (periodic). The periodicity depens on whether	the disksInCell parameter is given.

	*O.switchScene()* magic is used to have clean simulation for TriaxialTest without deleting the original simulation.
	This function therefore should never run in parallel with some code accessing your simulation.

	:param predicate: solid-defining predicate for which we generate packing
	:param disksInCell: if given, the packing will be periodic, with given number of disks in the periodic cell.
	:param radius: mean radius of disks
	:param rRelFuzz: relative fuzz of the radius -- e.g. radius=10, rRelFuzz=.2, then disks will have radii 10 ± (10*.2)), with an uniform distribution.
		0 by default, meaning all disks will have exactly the same radius.
	:param cropLayers: (aperiodic only) how many layers of disks will be added to the computed dimension of the box so that there no
		(or not so much, at least) boundary effects at the boundaries of the predicate.
	:param dim: dimension of the packing, to override dimensions of the predicate (if it is infinite, for instance)
	:param memoizeDb: name of sqlite database (existent or nonexistent) to find an already generated packing or to store
		the packing that will be generated, if not found (the technique of caching results of expensive computations
		is known as memoization). Fuzzy matching is used to select suitable candidate -- packing will be scaled, rRelFuzz
		and dimensions compared. Packing that are too small are dictarded. From the remaining candidate, the one with the
		least number disks will be loaded and returned.
	:param useOBB: effective only if a inGtsSurface predicate is given. If true (not default), oriented bounding box will be
		computed first; it can reduce substantially number of disks for the triaxial compression (like 10× depending on
		how much asymmetric the body is), see examples/gts-horse/gts-random-pack-obb.py
	:param memoDbg: show packings that are considered and reasons why they are rejected/accepted
	:param returnDiskPack: see the corresponding argument in :yref:`sudodem.pack.filterDiskPack`

	:return: DiskPack object with disks, filtered by the predicate.
	"""
	import sqlite3, os.path, cPickle, time, sys, _packPredicates, numpy
	from math import pi
	wantPeri=(disksInCell>0)
	if 'inGtsSurface' in dir(_packPredicates) and type(predicate)==inGtsSurface and useOBB:
		center,dim,orientation=gtsSurfaceBestFitOBB(predicate.surf)
		print "Best-fit oriented-bounding-box computed for GTS surface, orientation is",orientation
		dim*=2 # gtsSurfaceBestFitOBB returns halfSize
	else:
		if not dim: dim=predicate.dim()
		if max(dim)==float('inf'): raise RuntimeError("Infinite predicate and no dimension of packing requested.")
		center=predicate.center()
		orientation=None
	if not wantPeri: fullDim=tuple([dim[i]+4*cropLayers*radius for i in 0,1,2])
	else:
		# compute cell dimensions now, as they will be compared to ones stored in the db
		# they have to be adjusted to not make the cell to small WRT particle radius
		fullDim=dim
		cloudPorosity=0.25 # assume this number for the initial cloud (can be underestimated)
		beta,gamma=fullDim[1]/fullDim[0],fullDim[2]/fullDim[0] # ratios β=y₀/x₀, γ=z₀/x₀
		N100=disksInCell/cloudPorosity # number of disks for cell being filled by disks without porosity
		x1=radius*(1/(beta*gamma)*N100*(4/3.)*pi)**(1/3.)
		y1,z1=beta*x1,gamma*x1; vol0=x1*y1*z1
		maxR=radius*(1+rRelFuzz)
		x1=max(x1,8*maxR); y1=max(y1,8*maxR); z1=max(z1,8*maxR); vol1=x1*y1*z1
		N100*=vol1/vol0 # volume might have been increased, increase number of disks to keep porosity the same
		sp=_getMemoizedPacking(memoizeDb,radius,rRelFuzz,x1,y1,z1,fullDim,wantPeri,fillPeriodic=True,disksInCell=disksInCell,memoDbg=False)
		if sp:
			if orientation:
				sp.cellSize=(0,0,0) # resetting cellSize avoids warning when rotating
				sp.rotate(*orientation.toAxisAngle())
			return filterDiskPack(predicate,sp,material=material,returnDiskPack=returnDiskPack)
		else: print "No suitable packing in database found, running",'PERIODIC compression' if wantPeri else 'triaxial'
		sys.stdout.flush()
	O.switchScene(); O.resetThisScene() ### !!
	if wantPeri:
		# x1,y1,z1 already computed above
		sp=DiskPack()
		O.periodic=True
		#O.cell.refSize=(x1,y1,z1)
		O.cell.setBox((x1,y1,z1))
		#print cloudPorosity,beta,gamma,N100,x1,y1,z1,O.cell.refSize
		#print x1,y1,z1,radius,rRelFuzz
		O.materials.append(FrictMat(young=3e10,density=2400))
		num=sp.makeCloud(Vector3().Zero,O.cell.refSize,radius,rRelFuzz,disksInCell,True)
		O.engines=[ForceResetter(),InsertionSortCollider([Bo1_Disk_Aabb()],verletDist=.05*radius),InteractionLoop([Ig2_Disk_Disk_ScGeom()],[Ip2_FrictMat_FrictMat_FrictPhys()],[Law2_ScGeom_FrictPhys_CundallStrack()]),PeriIsoCompressor(charLen=2*radius,stresses=[-100e9,-1e8],maxUnbalanced=1e-2,doneHook='O.pause();',globalUpdateInt=5,keepProportions=True),NewtonIntegrator(damping=.6)]
		O.materials.append(FrictMat(young=30e9,frictionAngle=.5,poisson=.3,density=1e3))
		for s in sp: O.bodies.append(utils.disk(s[0],s[1]))
		O.dt=utils.PWaveTimeStep()
		O.run(); O.wait()
		sp=DiskPack(); sp.fromSimulation()
		#print 'Resulting cellSize',sp.cellSize,'proportions',sp.cellSize[1]/sp.cellSize[0],sp.cellSize[2]/sp.cellSize[0]
		# repetition to the required cell size will be done below, after memoizing the result
	else:
		assumedFinalDensity=0.6
		V=(4.0/3.0)*pi*radius**3.0; N=assumedFinalDensity*fullDim[0]*fullDim[1]*fullDim[2]/V;
		TriaxialTest(
			numberOfGrains=int(N),radiusMean=radius,radiusStdDev=rRelFuzz,
			# upperCorner is just size ratio, if radiusMean is specified
			upperCorner=fullDim,
			## no need to touch any the following
			noFiles=True,lowerCorner=[0,0,0],sigmaIsoCompaction=1e7,sigmaLateralConfinement=1e5,compactionFrictionDeg=1,StabilityCriterion=.05,strainRate=.2,thickness=-1,maxWallVelocity=.1,wallOversizeFactor=1.5,autoUnload=True,autoCompressionActivation=False).load()
		while ( numpy.isnan(utils.unbalancedForce()) or utils.unbalancedForce()>0.005 ) :
			O.run(100,True)
		sp=DiskPack(); sp.fromSimulation()
	O.switchScene() ### !!
	_memoizePacking(memoizeDb,sp,radius,rRelFuzz,wantPeri,fullDim)
	if wantPeri: sp.cellFill(Vector3(fullDim[0],fullDim[1],fullDim[2]))
	if orientation:
		sp.cellSize=(0,0,0); # reset periodicity to avoid warning when rotating periodic packing
		sp.rotate(*orientation.toAxisAngle())
	return filterDiskPack(predicate,sp,material=material,color=color,returnDiskPack=returnDiskPack)

def randomPeriPack(radius,initSize,rRelFuzz=0.0,memoizeDb=None,noPrint=False):
	"""Generate periodic dense packing.

	A cell of initSize is stuffed with as many disks as possible, then we run periodic compression with PeriIsoCompressor, just like with
	randomDensePack.

	:param radius: mean disk radius
	:param rRelFuzz: relative fuzz of disk radius (equal distribution); see the same param for randomDensePack.
	:param initSize: initial size of the periodic cell.

	:return: DiskPack object, which also contains periodicity information.
	"""
	from math import pi
	sp=_getMemoizedPacking(memoizeDb,radius,rRelFuzz,initSize[0],initSize[1],initSize[2],fullDim=Vector3(0,0,0),wantPeri=True,fillPeriodic=False,disksInCell=-1,memoDbg=True,noPrint=noPrint)
	if sp: return sp
	O.switchScene(); O.resetThisScene()
	sp=DiskPack()
	O.periodic=True
	#O.cell.refSize=initSize
	O.cell.setBox(initSize)
	sp.makeCloud(Vector3().Zero,O.cell.refSize,radius,rRelFuzz,-1,True)
	O.engines=[ForceResetter(),InsertionSortCollider([Bo1_Disk_Aabb()],verletDist=.05*radius),InteractionLoop([Ig2_Disk_Disk_ScGeom()],[Ip2_FrictMat_FrictMat_FrictPhys()],[Law2_ScGeom_FrictPhys_CundallStrack()]),PeriIsoCompressor(charLen=2*radius,stresses=[-100e9,-1e8],maxUnbalanced=1e-2,doneHook='O.pause();',globalUpdateInt=20,keepProportions=True),NewtonIntegrator(damping=.8)]
	O.materials.append(FrictMat(young=30e9,frictionAngle=.1,poisson=.3,density=1e3))
	for s in sp: O.bodies.append(utils.disk(s[0],s[1]))
	O.dt=utils.PWaveTimeStep()
	O.timingEnabled=True
	O.run(); O.wait()
	ret=DiskPack()
	ret.fromSimulation()
	_memoizePacking(memoizeDb,ret,radius,rRelFuzz,wantPeri=True,fullDim=Vector3(0,0,0),noPrint=noPrint) # fullDim unused
	O.switchScene()
	return ret

def hexaNet( radius, cornerCoord=[0,0,0], xLength=1., yLength=0.5, mos=0.08, a=0.04, b=0.04, startAtCorner=True, isSymmetric=False, **kw ):
	"""Definition of the particles for a hexagonal wire net in the x-y-plane for the WireMatPM.

	:param radius: radius of the particle
	:param cornerCoord: coordinates of the lower left corner of the net
	:param xLenght: net length in x-direction
	:param yLenght: net length in y-direction
	:param mos: mesh opening size (horizontal distance between the double twists)
	:param a: length of double-twist
	:param b: height of single wire section
	:param startAtCorner: if true the generation starts with a double-twist at the lower left corner
	:param isSymmetric: defines if the net is symmetric with respect to the y-axis

	:return: set of disks which defines the net (net) and exact dimensions of the net (lx,ly).

	note::
	This packing works for the WireMatPM only. The particles at the corner are always generated first. For examples on how to use this packing see examples/WireMatPM. In order to create the proper interactions for the net the interaction radius has to be adapted in the simulation.

	"""
	# check input dimension
	if(xLength<mos): raise ValueError("xLength must be greater than mos!");
	if(yLength<2*a+b): raise ValueError("yLength must be greater than 2*a+b!");
	xstart = cornerCoord[0]
	ystart = cornerCoord[1]
	z = cornerCoord[2]
	ab = a+b
	# number of double twisted sections in y-direction and real length ly
	ny = int( (yLength-a)/ab ) + 1
	ly = ny*a+(ny-1)*b
	jump=0
	# number of sections in x-direction and real length lx
	if isSymmetric:
		nx = int( xLength/mos ) + 1
		lx = (nx-1)*mos
		if not startAtCorner:
			nx+=-1
	else:
		nx = int( (xLength-0.5*mos)/mos ) + 1
		lx = (nx-1)*mos+0.5*mos
	net = []
	# generate corner particles
	if startAtCorner:
		if (ny%2==0): # if ny even no symmetry in y-direction
			net+=[utils.disk((xstart,ystart+ly,z),radius=radius,**kw)] # upper left corner
			if isSymmetric:
				net+=[utils.disk((xstart+lx,ystart+ly,z),radius=radius,**kw)] # upper right corner
			else:
				net+=[utils.disk((xstart+lx,ystart,z),radius=radius,**kw)] # lower right corner
		else: # if ny odd symmetry in y-direction
			if not isSymmetric:
				net+=[utils.disk((xstart+lx,ystart,z),radius=radius,**kw)] # lower right corner
				net+=[utils.disk((xstart+lx,ystart+ly,z),radius=radius,**kw)] # upper right corner
		jump=1
	else: # do not start at corner
		if (ny%2==0): # if ny even no symmetry in y-direction
			net+=[utils.disk((xstart,ystart,z),radius=radius,**kw)] # lower left corner
			if isSymmetric:
				net+=[utils.disk((xstart+lx,ystart,z),radius=radius,**kw)] # lower right corner
			else:
				net+=[utils.disk((xstart+lx,ystart+ly,z),radius=radius,**kw)] # upper right corner
		else: # if ny odd symmetry in y-direction
			net+=[utils.disk((xstart,ystart,z),radius=radius,**kw)] # lower left corner
			net+=[utils.disk((xstart,ystart+ly,z),radius=radius,**kw)] # upper left corner
			if isSymmetric:
				net+=[utils.disk((xstart+lx,ystart,z),radius=radius,**kw)] # lower right corner
				net+=[utils.disk((xstart+lx,ystart+ly,z),radius=radius,**kw)] # upper right corner
		xstart+=0.5*mos
	# generate other particles
	if isSymmetric:
		for i in range(ny):
			y = ystart + i*ab
			for j in range(nx):
				x = xstart + j*mos
				# add two particles of one vertical section (double-twist)
				net+=[utils.disk((x,y,z),radius=radius,**kw)]
				net+=[utils.disk((x,y+a,z),radius=radius,**kw)]
			# set values for next section
			xstart = xstart - 0.5*mos*pow(-1,i+jump)
			nx = int(nx + 1*pow(-1,i+jump))
	else:
		for i in range(ny):
			y = ystart + i*ab
			for j in range(nx):
				x = xstart + j*mos
				# add two particles of one vertical section (double-twist)
				net+=[utils.disk((x,y,z),radius=radius,**kw)]
				net+=[utils.disk((x,y+a,z),radius=radius,**kw)]
			# set values for next section
			xstart = xstart - 0.5*mos*pow(-1,i+jump)
	return [net,lx,ly]
