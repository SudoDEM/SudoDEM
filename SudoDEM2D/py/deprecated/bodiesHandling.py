# encoding: utf-8
"""
Miscellaneous functions, which are useful for handling bodies.
"""

from sudodem.wrapper import *
import utils,math,numpy

try:
	from minieigen import *
except ImportError:
	from miniEigen import *

#disksPackDimensions==================================================
def disksPackDimensions(idDisks=[],mask=-1):
	"""The function accepts the list of disks id's or list of bodies and calculates max and min dimensions, geometrical center.

	:param list idDisks: list of disks
	:param int mask: :yref:`Body.mask` for the checked bodies

	:return: dictionary with keys ``min`` (minimal dimension, Vector3), ``max`` (maximal dimension, Vector3), ``minId`` (minimal dimension disk Id, Vector3), ``maxId`` (maximal dimension disk Id, Vector3), ``center`` (central point of bounding box, Vector3), ``extends`` (sizes of bounding box, Vector3), ``volume`` (volume of disks, Real), ``mass`` (mass of disks, Real), ``number`` (number of disks, int),

	"""
	idDisksIter=[]

	if (len(idDisks)<1):
		#check mask
		ifSpherMask=[]
		if (mask>-1):   #The case, when only the mask was given, without list of ids
			for i in O.bodies:
				if ((i.mask&mask)<>0):
					ifSpherMask.append(i.id)
			if (len(ifSpherMask)<2):
				raise RuntimeWarning("Not enough bodies to analyze with given mask")
			else:
				idDisksIter=ifSpherMask
		else:
			raise RuntimeWarning("Only a list of particles with length > 1 can be analyzed")
	else:
		idDisksIter=idDisks


	minVal = Vector3.Zero
	maxVal = Vector3.Zero

	minId = Vector3.Zero
	maxId = Vector3.Zero

	counter = 0
	volume = 0.0
	mass = 0.0


	for i in idDisksIter:
		if (type(i).__name__=='int'):
			b = O.bodies[i]			#We have received a list of ID's
		elif (type(i).__name__=='Body'):
			b = i								#We have recevied a list of bodies
		else:
			raise TypeError("Unknow type of data, should be list of int's or bodies's")

		if (b):
			diskPosition=b.state.pos	#skip non-existent disks

			try:
				diskRadius=b.shape.radius	#skip non-disks
			except AttributeError: continue

			if (mask>-1) and ((mask&b.mask)==0): continue			#skip bodies with wrong mask


			diskRadiusVec3 = Vector3(diskRadius,diskRadius,diskRadius)

			diskMax = diskPosition + diskRadiusVec3
			diskMin = diskPosition - diskRadiusVec3

			for dim in range(0,3):
				if ((diskMax[dim]>maxVal[dim]) or (counter==0)):
					maxVal[dim]=diskMax[dim]
					maxId[dim] = b.id
				if ((diskMin[dim]<minVal[dim]) or (counter==0)):
					minVal[dim]=diskMin[dim]
					minId[dim] = b.id
			volume += 4.0/3.0*math.pi*diskRadius*diskRadius*diskRadius
			mass += b.state.mass
			counter += 1

	center = (maxVal-minVal)/2.0+minVal
	extends = maxVal-minVal

	dimensions = {'max':maxVal,'min':minVal,'maxId':maxId,'minId':minId,'center':center,
		'extends':extends, 'volume':volume, 'mass':mass, 'number':counter}
	return dimensions

#facetsDimensions==================================================
def facetsDimensions(idFacets=[],mask=-1):
	"""The function accepts the list of facet id's or list of facets and calculates max and min dimensions, geometrical center.

	:param list idFacets: list of disks
	:param int mask: :yref:`Body.mask` for the checked bodies

	:return: dictionary with keys ``min`` (minimal dimension, Vector3), ``max`` (maximal dimension, Vector3), ``minId`` (minimal dimension facet Id, Vector3), ``maxId`` (maximal dimension facet Id, Vector3), ``center`` (central point of bounding box, Vector3), ``extends`` (sizes of bounding box, Vector3), ``number`` (number of facets, int),

	"""
	idFacetsIter=[]

	if (len(idFacets)<1):
		#check mask
		ifFacetMask=[]
		if (mask>-1):   #The case, when only the mask was given, without list of ids
			for i in O.bodies:
				if ((i.mask&mask)<>0):
					ifFacetMask.append(i.id)
			if (len(ifFacetMask)<2):
				raise RuntimeWarning("Not enough bodies to analyze with given mask")
			else:
				idFacetsIter=ifFacetMask
		else:
			raise RuntimeWarning("Only a list of particles with length > 1 can be analyzed")
	else:
		idFacetsIter=idFacets


	minVal = Vector3.Zero
	maxVal = Vector3.Zero

	minId = Vector3.Zero
	maxId = Vector3.Zero

	counter = 0


	for i in idFacetsIter:
		if (type(i).__name__=='int'):
			b = O.bodies[i]			#We have received a list of ID's
		elif (type(i).__name__=='Body'):
			b = i								#We have recevied a list of bodies
		else:
			raise TypeError("Unknow type of data, should be list of int's or bodies's")

		if (b):
			p = b.state.pos
			o = b.state.ori
			s = b.shape
			pt1 = p + o*s.vertices[0]
			pt2 = p + o*s.vertices[1]
			pt3 = p + o*s.vertices[2]

			if (mask>-1) and ((mask&b.mask)==0): continue			#skip bodies with wrong mask

			facetMax = Vector3(max(pt1[0], pt2[0], pt3[0]), max(pt1[1], pt2[1], pt3[1]), max(pt1[2], pt2[2], pt3[2]))
			facetMin = Vector3(min(pt1[0], pt2[0], pt3[0]), min(pt1[1], pt2[1], pt3[1]), min(pt1[2], pt2[2], pt3[2]))

			for dim in range(0,3):
				if ((facetMax[dim]>maxVal[dim]) or (counter==0)):
					maxVal[dim]=facetMax[dim]
					maxId[dim] = b.id
				if ((facetMin[dim]<minVal[dim]) or (counter==0)):
					minVal[dim]=facetMin[dim]
					minId[dim] = b.id
			counter += 1

	center = (maxVal-minVal)/2.0+minVal
	extends = maxVal-minVal

	dimensions = {'max':maxVal,'min':minVal,'maxId':maxId,'minId':minId,'center':center,
		'extends':extends, 'number':counter}
	return dimensions

#disksPackDimensions==================================================
def disksModify(idDisks=[],mask=-1,shift=Vector3.Zero,scale=1.0,orientation=Quaternion((0,1,0),0.0),copy=False):
	"""The function accepts the list of disks id's or list of bodies and modifies them: rotating, scaling, shifting.
	if copy=True copies bodies and modifies them.
	Also the mask can be given. If idDisks not empty, the function affects only bodies, where the mask passes.
	If idDisks is empty, the function search for bodies, where the mask passes.

	:param Vector3 shift: Vector3(X,Y,Z) parameter moves disks.
	:param float scale: factor scales given disks.
	:param Quaternion orientation: orientation of disks
	:param int mask: :yref:`Body.mask` for the checked bodies
	:returns: list of bodies if copy=True, and Boolean value if copy=False
	"""

	idDisksIter=[]

	if (len(idDisks)==0):
		#check mask
		ifSpherMask=[]
		if (mask>-1):   #The case, when only the mask was given, without list of ids
			for i in O.bodies:
				if ((i.mask&mask)<>0):
					ifSpherMask.append(i.id)
			if (len(ifSpherMask)==0):
				raise RuntimeWarning("No bodies to modify with given mask")
			else:
				idDisksIter=ifSpherMask
		else:
			raise RuntimeWarning("No bodies to modify")
	else:
		idDisksIter=idDisks

	dims = disksPackDimensions(idDisksIter)

	ret=[]
	for i in idDisksIter:
		if (type(i).__name__=='int'):
			b = O.bodies[i]			#We have received a list of ID's
		elif (type(i).__name__=='Body'):
			b = i								#We have recevied a list of bodies
		else:
			raise TypeError("Unknown type of data, should be list of int's or bodies")

		try:
			diskRadius=b.shape.radius	#skip non-disks
		except AttributeError: continue

		if (mask>-1) and ((mask&b.mask)==0): continue			#skip bodies with wrong mask

		if (copy): b=diskDuplicate(b)

		b.state.pos=orientation*(b.state.pos-dims['center'])+dims['center']
		b.shape.radius*=scale
		b.state.pos=(b.state.pos-dims['center'])*scale + dims['center']

		b.state.pos+=shift

		if (copy): ret.append(b)

	if (copy):
		return ret
	else:
		return True

#disksDublicate=======================================================
def diskDuplicate(idDisk):
	"""The functions makes a copy of disk"""

	i=idDisk
	if (type(i).__name__=='int'):
		b = O.bodies[i]			#We have received a list of ID's
	elif (type(i).__name__=='Body'):
		b = i								#We have recevied a list of bodies
	else:
		raise TypeError("Unknown type of data, should be list of int's or bodies")

	try:
		diskRadius=b.shape.radius	#skip non-disks
	except AttributeError:
		return False

	addedBody = utils.disk(center=b.state.pos,radius=b.shape.radius,fixed=not(b.dynamic),wire=b.shape.wire,color=b.shape.color,highlight=b.shape.highlight,material=b.material,mask=b.mask)

	return addedBody


