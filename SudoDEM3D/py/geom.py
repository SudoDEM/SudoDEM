# encoding: utf-8
"""
Creates geometry objects from facets.
"""

from sudodem.wrapper import *
from sudodem import utils
import math,numpy

try:
	from minieigen import *
except ImportError:
	from miniEigen import *

#facetBox===============================================================
def facetBox(center,extents,orientation=Quaternion((0,1,0),0.0),wallMask=63,**kw):
	"""
	Create arbitrarily-aligned box composed of facets, with given center, extents and orientation.
	If any of the box dimensions is zero, corresponding facets will not be created. The facets are oriented outwards from the box.

	:param Vector3 center: center of the box
	:param Vector3 extents: lengths of the box sides
	:param Quaternion orientation: orientation of the box
	:param bitmask wallMask: determines which walls will be created, in the order -x (1), +x (2), -y (4), +y (8), -z (16), +z (32). The numbers are ANDed; the default 63 means to create all walls
	:param \*\*kw: (unused keyword arguments) passed to :yref:`sudodem.utils.facet`
	:returns: list of facets forming the box
	"""

	return facetParallelepiped(center=center, extents=extents, height=extents[2], orientation=orientation, wallMask=wallMask, **kw)

#facetParallelepiped===============================================================
def facetParallelepiped(center,extents,height,orientation=Quaternion((0,1,0),0.0),wallMask=63,**kw):
	"""
	Create arbitrarily-aligned Parallelepiped composed of facets, with given center, extents, height  and orientation.
	If any of the parallelepiped dimensions is zero, corresponding facets will not be created. The facets are oriented outwards from the parallelepiped.

	:param Vector3 center: center of the parallelepiped
	:param Vector3 extents: lengths of the parallelepiped sides
	:param Real height: height of the parallelepiped (along axis z)
	:param Quaternion orientation: orientation of the parallelepiped
	:param bitmask wallMask: determines which walls will be created, in the order -x (1), +x (2), -y (4), +y (8), -z (16), +z (32). The numbers are ANDed; the default 63 means to create all walls
	:param \*\*kw: (unused keyword arguments) passed to :yref:`sudodem.utils.facet`
	:returns: list of facets forming the parallelepiped
	"""

	if (height<0): raise RuntimeError("The height should have the positive value");
	if (height>extents[2]): raise RuntimeError("The height should be smaller or equal as extents[2]");

	#Defense from zero dimensions
	if (wallMask>63):
		print("wallMask must be 63 or less")
		wallMask=63
	if (extents[0]==0):
		wallMask=1
	elif (extents[1]==0):
		wallMask=4
	elif (extents[2]==0 or height==0):
		wallMask=16
	if (((extents[0]==0) and (extents[1]==0)) or ((extents[0]==0) and (extents[2]==0)) or ((extents[1]==0) and (extents[2]==0))):
		raise RuntimeError("Please, specify at least 2 none-zero dimensions in extents!");
	# ___________________________

	#inclination angle
	beta = 0; dx = 0
	if (height>0):
		beta = math.asin(height/extents[2])
		dx = math.cos(beta)*extents[2]

	mn,mx=[-extents[i] for i in [0,1,2]],[extents[i] for i in [0,1,2]]
	def doWall(a,b,c,d):
		return [utils.facet((a,b,c),**kw),utils.facet((a,c,d),**kw)]
	ret=[]

	mn[2] = -height
	mx[2] = +height

	A=orientation*Vector3(mn[0],mn[1],mn[2])+center
	B=orientation*Vector3(mx[0],mn[1],mn[2])+center
	C=orientation*Vector3(mx[0],mx[1],mn[2])+center
	D=orientation*Vector3(mn[0],mx[1],mn[2])+center
	E=orientation*Vector3(mn[0]+dx,mn[1],mx[2])+center
	F=orientation*Vector3(mx[0]+dx,mn[1],mx[2])+center
	G=orientation*Vector3(mx[0]+dx,mx[1],mx[2])+center
	H=orientation*Vector3(mn[0]+dx,mx[1],mx[2])+center
	if wallMask&1:  ret+=doWall(A,D,H,E)
	if wallMask&2:  ret+=doWall(B,F,G,C)
	if wallMask&4:  ret+=doWall(A,E,F,B)
	if wallMask&8:  ret+=doWall(D,C,G,H)
	if wallMask&16: ret+=doWall(A,B,C,D)
	if wallMask&32: ret+=doWall(E,H,G,F)
	return ret

#facetCylinder==========================================================
def facetCylinder(center,radius,height,orientation=Quaternion((0,1,0),0.0),
	segmentsNumber=10,wallMask=7,angleRange=None,closeGap=False,
	radiusTopInner=-1, radiusBottomInner=-1,
	**kw):
	"""
	Create arbitrarily-aligned cylinder composed of facets, with given center, radius, height and orientation.
	Return List of facets forming the cylinder;

	:param Vector3 center: center of the created cylinder
	:param float radius:  cylinder radius
	:param float height: cylinder height
	:param float radiusTopInner: inner radius of cylinders top, -1 by default
	:param float radiusBottomInner: inner radius of cylinders bottom, -1 by default
	:param Quaternion orientation: orientation of the cylinder; the reference orientation has axis along the $+x$ axis.
	:param int segmentsNumber: number of edges on the cylinder surface (>=5)
	:param bitmask wallMask: determines which walls will be created, in the order up (1), down (2), side (4). The numbers are ANDed; the default 7 means to create all walls
	:param (θmin,Θmax) angleRange: allows one to create only part of bunker by specifying range of angles; if ``None``, (0,2*pi) is assumed.
	:param bool closeGap: close range skipped in angleRange with triangular facets at cylinder bases.
	:param \*\*kw: (unused keyword arguments) passed to utils.facet;
	"""
	# check zero dimentions
	if (radius<=0): raise RuntimeError("The radius should have the positive value");
	if (height<=0): wallMask = 1;

	return facetCylinderConeGenerator(center=center,radiusTop=radius,height=height,
		orientation=orientation,segmentsNumber=segmentsNumber,wallMask=wallMask,
		angleRange=angleRange,closeGap=closeGap,
		radiusTopInner=radiusTopInner, radiusBottomInner=radiusBottomInner,
		**kw)

#facetSphere==========================================================
def facetSphere(center,radius,thetaResolution=8,phiResolution=8,returnElementMap=False,**kw):
	"""
	Create arbitrarily-aligned sphere composed of facets, with given center, radius and orientation.
	Return List of facets forming the sphere. Parameters inspired by ParaView sphere glyph

	:param Vector3 center: center of the created sphere
	:param float radius: sphere radius
	:param int thetaResolution: number of facets around "equator"
	:param int phiResolution: number of facets between "poles" + 1
	:param bool returnElementMap: returns also tuple of nodes ((x1,y1,z1),(x2,y2,z2),...) and elements ((id01,id02,id03),(id11,id12,id13),...) if true, only facets otherwise
	:param \*\*kw: (unused keyword arguments) passed to utils.facet;
	"""
	# check zero dimentions
	if (radius<=0):          raise RuntimeError("The radius should have the positive value");
	if (thetaResolution<3): raise RuntimeError("thetaResolution must be > 3");
	if (phiResolution<3):   raise RuntimeError("phiResolution must be > 3");

	r,c0,c1,c2 = radius,center[0],center[1],center[2]
	nodes = [Vector3(c0,c1,c2+radius)]
	phis   = numpy.linspace(math.pi/(phiResolution-1),math.pi,phiResolution-2,endpoint=False)
	thetas = numpy.linspace(0,2*math.pi,thetaResolution,endpoint=False)
	nodes.extend((Vector3(c0+r*math.cos(theta)*math.sin(phi),c1+r*math.sin(theta)*math.sin(phi),c2+r*math.cos(phi)) for phi in phis for theta in thetas))
	nodes.append(Vector3(c0,c1,c2-radius))
	n = len(nodes)-1

	elements = [(0,i+1,i+2) for i in xrange(thetaResolution-1)]
	elements.append((0,1,thetaResolution))
	for j in xrange(0,phiResolution-3):
		k = j*thetaResolution + 1
		elements.extend((k+i,k+i+1,k+i+thetaResolution) for i in xrange(thetaResolution-1))
		elements.append((k,k+thetaResolution-1,k+2*thetaResolution-1))
		elements.extend((k+i+thetaResolution,k+i+1+thetaResolution,k+i+1) for i in xrange(thetaResolution-1))
		elements.append((k+2*thetaResolution-1,k+thetaResolution,k))
	elements.extend((n,n-i-1,n-i-2) for i in xrange(thetaResolution-1))
	elements.append((n,n-1,n-thetaResolution))

	facets = [utils.facet(tuple(nodes[node] for node in elem),**kw) for elem in elements]
	if returnElementMap:
		return facets,nodes,elements
	return facets


#facetCone==============================================================
def facetCone(center,radiusTop,radiusBottom,height,orientation=Quaternion((0,1,0),0.0),
	segmentsNumber=10,wallMask=7,angleRange=None,closeGap=False,
	radiusTopInner=-1, radiusBottomInner=-1,
	**kw):
	"""
	Create arbitrarily-aligned cone composed of facets, with given center, radius, height and orientation.
	Return List of facets forming the cone;

	:param Vector3 center: center of the created cylinder
	:param float radiusTop:  cone top radius
	:param float radiusBottom:  cone bottom radius
	:param float radiusTopInner: inner radius of cones top, -1 by default
	:param float radiusBottomInner: inner radius of cones bottom, -1 by default
	:param float height: cone height
	:param Quaternion orientation: orientation of the cone; the reference orientation has axis along the $+x$ axis.
	:param int segmentsNumber: number of edges on the cone surface (>=5)
	:param bitmask wallMask: determines which walls will be created, in the order up (1), down (2), side (4). The numbers are ANDed; the default 7 means to create all walls
	:param (θmin,Θmax) angleRange: allows one to create only part of cone by specifying range of angles; if ``None``, (0,2*pi) is assumed.
	:param bool closeGap: close range skipped in angleRange with triangular facets at cylinder bases.
	:param \*\*kw: (unused keyword arguments) passed to utils.facet;
	"""
	# check zero dimentions
	if ((radiusBottom<=0) and (radiusTop<=0)): raise RuntimeError("The radiusBottom or radiusTop should have the positive value");

	return facetCylinderConeGenerator(center=center,radiusTop=radiusTop,
		radiusBottom=radiusBottom,height=height,orientation=orientation,segmentsNumber=segmentsNumber,
		wallMask=wallMask,angleRange=angleRange,closeGap=closeGap,
		radiusTopInner=radiusTopInner, radiusBottomInner=radiusBottomInner,
		**kw)

#facetPolygon===========================================================
def facetPolygon(center,radiusOuter,orientation=Quaternion((0,1,0),0.0),segmentsNumber=10,angleRange=None,radiusInner=0,**kw):
	"""
	Create arbitrarily-aligned polygon composed of facets, with given center, radius (outer and inner) and orientation.
	Return List of facets forming the polygon;

	:param Vector3 center: center of the created cylinder
	:param float radiusOuter:  outer radius
	:param float radiusInner: inner height (can be 0)
	:param Quaternion orientation: orientation of the polygon; the reference orientation has axis along the $+x$ axis.
	:param int segmentsNumber: number of edges on the polygon surface (>=3)
	:param (θmin,Θmax) angleRange: allows one to create only part of polygon by specifying range of angles; if ``None``, (0,2*pi) is assumed.
	:param \*\*kw: (unused keyword arguments) passed to utils.facet;
	"""
	# check zero dimentions
	if (abs(angleRange[1]-angleRange[0])>2.0*math.pi): raise RuntimeError("The |angleRange| cannot be larger 2.0*math.pi");

	return facetPolygonHelixGenerator(center=center,radiusOuter=radiusOuter,orientation=orientation,segmentsNumber=segmentsNumber,angleRange=angleRange,radiusInner=radiusInner,**kw)

#facetHelix===========================================================
def facetHelix(center,radiusOuter,pitch,orientation=Quaternion((0,1,0),0.0),segmentsNumber=10,angleRange=None,radiusInner=0,**kw):
	"""
	Create arbitrarily-aligned helix composed of facets, with given center, radius (outer and inner), pitch and orientation.
	Return List of facets forming the helix;

	:param Vector3 center: center of the created cylinder
	:param float radiusOuter:  outer radius
	:param float radiusInner: inner height (can be 0)
	:param Quaternion orientation: orientation of the helix; the reference orientation has axis along the $+x$ axis.
	:param int segmentsNumber: number of edges on the helix surface (>=3)
	:param (θmin,Θmax) angleRange: range of angles; if ``None``, (0,2*pi) is assumed.
	:param \*\*kw: (unused keyword arguments) passed to utils.facet;
	"""

	# check zero dimentions
	if (pitch<=0): raise RuntimeError("The pitch should have the positive value");
	return facetPolygonHelixGenerator(center=center,radiusOuter=radiusOuter,orientation=orientation,segmentsNumber=segmentsNumber,angleRange=angleRange,radiusInner=radiusInner,pitch=pitch,**kw)

#facetBunker============================================================
def facetBunker(center,dBunker,dOutput,hBunker,hOutput,hPipe=0.0,orientation=Quaternion((0,1,0),0.0),segmentsNumber=10,wallMask=4,angleRange=None,closeGap=False,**kw):
	"""
	Create arbitrarily-aligned bunker, composed of facets, with given center, radii, heights and orientation.
	Return List of facets forming the bunker;

	.. code-block:: none

		   dBunker
		______________
		|            |
		|            |
		|            | hBunker
		|            |
		|            |
		|            |
		|____________|
		\            /
		 \          /
		  \        /   hOutput
		   \      /
		    \____/
		    |    |
		    |____|     hPipe
		    dOutput

	:param Vector3 center: center of the created bunker
	:param float dBunker: bunker diameter, top
	:param float dOutput: bunker output diameter
	:param float hBunker: bunker height
	:param float hOutput: bunker output height
	:param float hPipe: bunker pipe height
	:param Quaternion orientation: orientation of the bunker; the reference orientation has axis along the $+x$ axis.
	:param int segmentsNumber: number of edges on the bunker surface (>=5)
	:param bitmask wallMask: determines which walls will be created, in the order up (1), down (2), side (4). The numbers are ANDed; the default 7 means to create all walls
	:param (θmin,Θmax) angleRange: allows one to create only part of bunker by specifying range of angles; if ``None``, (0,2*pi) is assumed.
	:param bool closeGap: close range skipped in angleRange with triangular facets at cylinder bases.
	:param \*\*kw: (unused keyword arguments) passed to utils.facet;
	"""
	# check zero dimentions
	if (dBunker<=0): raise RuntimeError("The diameter dBunker should have the positive value");
	if (dOutput<=0): raise RuntimeError("The diameter dOutput should have the positive value");
	if (hBunker<0): raise RuntimeError("The height hBunker should have the positive or or zero");
	if (hOutput<=0): raise RuntimeError("The height hOutput should have the positive value");
	if (hPipe<0): raise RuntimeError("The height hPipe should have the positive value or zero");

	ret=[]
	if ((hPipe>0) or (wallMask&2)):
		centerPipe = Vector3(0,0,hPipe/2.0)
		ret+=facetCylinder(center=centerPipe,radius=dOutput/2.0,height=hPipe,segmentsNumber=segmentsNumber,wallMask=wallMask&6,angleRange=angleRange,closeGap=closeGap,**kw)

	centerOutput = Vector3(0.0,0.0,hPipe+hOutput/2.0)
	ret+=facetCone(center=centerOutput,radiusTop=dBunker/2.0,radiusBottom=dOutput/2.0,height=hOutput,segmentsNumber=segmentsNumber,wallMask=wallMask&4,angleRange=angleRange,closeGap=closeGap,**kw)

	if (hBunker>0):
		centerBunker = Vector3(0.0,0.0,hPipe+hOutput+hBunker/2.0)
		ret+=facetCylinder(center=centerBunker,radius=dBunker/2.0,height=hBunker,segmentsNumber=segmentsNumber,wallMask=wallMask&5,angleRange=angleRange,closeGap=closeGap,**kw)

	for i in ret:
		i.state.pos=orientation*(i.state.pos)+Vector3(center)
		i.state.ori=orientation

	return ret

#facetPolygonHelixGenerator==================================================
def facetPolygonHelixGenerator(center,radiusOuter,pitch=0,orientation=Quaternion((0,1,0),0.0),segmentsNumber=10,angleRange=None,radiusInner=0,**kw):
	"""
	Please, do not use this function directly! Use geom.facetPloygon and geom.facetHelix instead.
	This is the base function for generating polygons and helixes from facets.
	"""
	# check zero dimentions
	if (segmentsNumber<3): raise RuntimeError("The segmentsNumber should be at least 3");
	if (radiusOuter<=0): raise RuntimeError("The radiusOuter should have the positive value");
	if (radiusInner<0): raise RuntimeError("The radiusInner should have the positive value or 0");
	if angleRange==None: angleRange=(0,2*math.pi)

	anglesInRad = numpy.linspace(angleRange[0], angleRange[1], segmentsNumber+1, endpoint=True)
	heightsInRad = numpy.linspace(0, pitch*(abs(angleRange[1]-angleRange[0])/(2.0*math.pi)), segmentsNumber+1, endpoint=True)

	POuter=[];
	PInner=[];
	PCenter=[];
	z=0;
	for i in anglesInRad:
		XOuter=radiusOuter*math.cos(i); YOuter=radiusOuter*math.sin(i);
		POuter.append(Vector3(XOuter,YOuter,heightsInRad[z]))
		PCenter.append(Vector3(0,0,heightsInRad[z]))
		if (radiusInner!=0):
			XInner=radiusInner*math.cos(i); YInner=radiusInner*math.sin(i);
			PInner.append(Vector3(XInner,YInner,heightsInRad[z]))
		z+=1

	for i in range(0,len(POuter)):
		POuter[i]=orientation*POuter[i]+center
		PCenter[i]=orientation*PCenter[i]+center
		if (radiusInner!=0):
			PInner[i]=orientation*PInner[i]+center

	ret=[]
	for i in range(1,len(POuter)):
		if (radiusInner==0):
			ret.append(utils.facet((PCenter[i],POuter[i],POuter[i-1]),**kw))
		else:
			ret.append(utils.facet((PInner[i-1],POuter[i-1],POuter[i]),**kw))
			ret.append(utils.facet((PInner[i],PInner[i-1],POuter[i]),**kw))

	return ret


#facetCylinderConeGenerator=============================================
def facetCylinderConeGenerator(center,radiusTop,height,orientation=Quaternion((0,1,0),0.0),
	segmentsNumber=10,wallMask=7,angleRange=None,closeGap=False,
	radiusBottom=-1,
	radiusTopInner=-1,
	radiusBottomInner=-1,
	**kw):
	"""
	Please, do not use this function directly! Use geom.facetCylinder and geom.facetCone instead.
	This is the base function for generating cylinders and cones from facets.
	:param float radiusTop:  top radius
	:param float radiusBottom:  bottom radius
	:param \*\*kw: (unused keyword arguments) passed to utils.facet;
	"""

	#For cylinders top and bottom radii are equal
	if (radiusBottom == -1):
		radiusBottom = radiusTop

	if ((radiusTopInner > 0 and radiusTopInner > radiusTop) or (radiusBottomInner > 0 and radiusBottomInner > radiusBottom)):
		raise RuntimeError("The internal radius cannot be larger than outer");
	# check zero dimentions
	if (segmentsNumber<3): raise RuntimeError("The segmentsNumber should be at least 3");
	if (height<0): raise RuntimeError("The height should have the positive value");
	if angleRange==None: angleRange=(0,2*math.pi)
	if (abs(angleRange[1]-angleRange[0])>2.0*math.pi): raise RuntimeError("The |angleRange| cannot be larger 2.0*math.pi");
	if (angleRange[1]<angleRange[0]): raise RuntimeError("angleRange[1] should be larger or equal angleRange[1]");

	if isinstance(angleRange,float):
		print(u'WARNING: geom.facetCylinder,angleRange should be (Θmin,Θmax), not just Θmax (one number), update your code.')
		angleRange=(0,angleRange)

	anglesInRad = numpy.linspace(angleRange[0], angleRange[1], segmentsNumber+1, endpoint=True)

	PTop=[]; PTop.append(Vector3(0,0,+height/2))
	PTopIn=[]; PTopIn.append(Vector3(0,0,+height/2))

	PBottom=[]; PBottom.append(Vector3(0,0,-height/2))
	PBottomIn=[]; PBottomIn.append(Vector3(0,0,-height/2))

	for i in anglesInRad:
		XTop=radiusTop*math.cos(i); YTop=radiusTop*math.sin(i);
		PTop.append(Vector3(XTop,YTop,+height/2))
		if (radiusTopInner > 0):
			XTopIn=radiusTopInner*math.cos(i); YTopIn=radiusTopInner*math.sin(i);
			PTopIn.append(Vector3(XTopIn,YTopIn,+height/2))

		XBottom=radiusBottom*math.cos(i); YBottom=radiusBottom*math.sin(i);
		PBottom.append(Vector3(XBottom,YBottom,-height/2))
		if (radiusBottomInner > 0):
			XBottomIn=radiusBottomInner*math.cos(i); YBottomIn=radiusBottomInner*math.sin(i);
			PBottomIn.append(Vector3(XBottomIn,YBottomIn,-height/2))

	for i in range(0,len(PTop)):
		PTop[i]=orientation*PTop[i]+center
		PBottom[i]=orientation*PBottom[i]+center
		if (len(PTopIn)>1):
			PTopIn[i]=orientation*PTopIn[i]+center
		if (len(PBottomIn)>1):
			PBottomIn[i]=orientation*PBottomIn[i]+center

	ret=[]
	for i in range(2,len(PTop)):
		if (wallMask&1)and(radiusTop!=0):
			if (len(PTopIn)>1):
				ret.append(utils.facet((PTop[i-1],PTopIn[i],PTopIn[i-1]),**kw))
				ret.append(utils.facet((PTop[i-1],PTop[i],PTopIn[i]),**kw))
			else:
				ret.append(utils.facet((PTop[0],PTop[i],PTop[i-1]),**kw))

		if (wallMask&2)and(radiusBottom!=0):
			if (len(PBottomIn)>1):
				ret.append(utils.facet((PBottom[i-1],PBottomIn[i],PBottomIn[i-1]),**kw))
				ret.append(utils.facet((PBottom[i-1],PBottom[i],PBottomIn[i]),**kw))
			else:
				ret.append(utils.facet((PBottom[0],PBottom[i-1],PBottom[i]),**kw))

		if wallMask&4:
			if (radiusBottom!=0):
				ret.append(utils.facet((PTop[i],PBottom[i],PBottom[i-1]),**kw))
			if (radiusTop!=0):
				ret.append(utils.facet((PBottom[i-1],PTop[i-1],PTop[i]),**kw))

	if (closeGap):
		if (wallMask&1)and(radiusTop!=0)and(abs(((angleRange[1]-angleRange[0])) > math.pi)):
			pts=[(radiusTop*math.cos(angleRange[i]),radiusTop*math.sin(angleRange[i])) for i in (0,1)]
			pp=[(pts[0][0],pts[0][1],+height/2.0), (pts[1][0],pts[1][1],+height/2.0), (0,0,+height/2.0)]
			pp=[orientation*p+center for p in pp]
			ret.append(utils.facet(pp,**kw))

		if (wallMask&2)and(radiusBottom!=0)and(abs(((angleRange[1]-angleRange[0])) > math.pi)):
			pts=[(radiusBottom*math.cos(angleRange[i]),radiusBottom*math.sin(angleRange[i])) for i in (0,1)]
			pp=[(0,0,-height/2.0), (pts[1][0],pts[1][1],-height/2.0), (pts[0][0],pts[0][1],-height/2.0)]
			pp=[orientation*p+center for p in pp]
			ret.append(utils.facet(pp,**kw))

		if (wallMask&4):
			ptsBottom=[(radiusBottom*math.cos(angleRange[i]),radiusBottom*math.sin(angleRange[i])) for i in (0,1)]
			ptsTop=[(radiusTop*math.cos(angleRange[i]),radiusTop*math.sin(angleRange[i])) for i in (0,1)]

			if (abs(((angleRange[1]-angleRange[0])) >= math.pi)):
				if (radiusBottom!=0)and(radiusTop!=0):	#Cylinder
					pp=[(ptsBottom[0][0],ptsBottom[0][1],-height/2.0),(ptsBottom[1][0],ptsBottom[1][1],-height/2.0),(ptsTop[0][0],ptsTop[0][1],height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))

					pp=[(ptsBottom[1][0],ptsBottom[1][1],-height/2.0), (ptsTop[1][0],ptsTop[1][1],height/2.0), (ptsTop[0][0],ptsTop[0][1],height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))
				elif (radiusBottom==0)and(radiusTop!=0):	#ConeTop
					pp=[(ptsTop[1][0],ptsTop[1][1],height/2.0), (ptsTop[0][0],ptsTop[0][1],height/2.0), (0,0,-height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))
				elif (radiusTop==0)and(radiusBottom!=0):	#ConeBottom
					pp=[(0,0,height/2.0),(ptsBottom[0][0],ptsBottom[0][1],-height/2.0),(ptsBottom[1][0],ptsBottom[1][1],-height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))
			else:
				if (radiusBottom!=0)and(radiusTop!=0):	#Cylinder
					pp=[(ptsBottom[0][0],ptsBottom[0][1],-height/2.0),(0,0,-height/2.0),(ptsTop[0][0],ptsTop[0][1],height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))

					pp=[(0,0,-height/2.0), (0,0,height/2.0), (ptsTop[0][0],ptsTop[0][1],height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))

					pp=[(0,0,-height/2.0),(ptsBottom[1][0],ptsBottom[1][1],-height/2.0),(0,0,height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))

					pp=[(ptsBottom[1][0],ptsBottom[1][1],-height/2.0), (ptsTop[1][0],ptsTop[1][1],height/2.0), (0,0,height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))
				elif (radiusBottom==0)and(radiusTop!=0):	#ConeTop
					pp=[(0,0,height/2.0), (ptsTop[0][0],ptsTop[0][1],height/2.0), (0,0,-height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))

					pp=[(ptsTop[1][0],ptsTop[1][1],height/2.0), (0,0,height/2.0), (0,0,-height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))
				elif (radiusTop==0)and(radiusBottom!=0):	#ConeBottom
					pp=[(0,0,height/2.0),(ptsBottom[0][0],ptsBottom[0][1],-height/2.0),(0,0,-height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))

					pp=[(0,0,height/2.0),(0,0,-height/2.0),(ptsBottom[1][0],ptsBottom[1][1],-height/2.0)]
					pp=[orientation*p+center for p in pp]
					ret.append(utils.facet(pp,**kw))
	return ret
