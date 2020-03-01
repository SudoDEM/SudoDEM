"""
Import geometry from various formats ('import' is python keyword, hence the name 'ymport').
"""

from sudodem.wrapper import *
from sudodem import utils

try:
	from minieigen import *
except ImportError:
	from miniEigen import *

def textExt(fileName,format='x_y_z_r',shift=Vector3.Zero,scale=1.0,**kw):
	"""Load sphere coordinates from file in specific format, returns a list of corresponding bodies; that may be inserted to the simulation with O.bodies.append().

	:param str filename: file name
	:param str format: the name of output format. Supported `x_y_z_r`(default), `x_y_z_r_matId`
	:param [float,float,float] shift: [X,Y,Z] parameter moves the specimen.
	:param float scale: factor scales the given data.
	:param \*\*kw: (unused keyword arguments) is passed to :yref:`sudodem.utils.sphere`
	:returns: list of spheres.

	Lines starting with # are skipped
	"""
	infile = open(fileName,"r")
	lines = infile.readlines()
	infile.close()
	ret=[]
	for line in lines:
		data = line.split()
		if (data[0] == "#format"):
			format=data[1]
			continue
		elif (data[0][0] == "#"): continue

		if (format=='x_y_z_r'):
			pos = Vector3(float(data[0]),float(data[1]),float(data[2]))
			ret.append(utils.sphere(shift+scale*pos,scale*float(data[3]),**kw))
		elif (format=='x_y_z_r_matId'):
			pos = Vector3(float(data[0]),float(data[1]),float(data[2]))
			ret.append(utils.sphere(shift+scale*pos,scale*float(data[3]),material=int(data[4]),**kw))

		elif (format=='id_x_y_z_r_matId'):
			pos = Vector3(float(data[1]),float(data[2]),float(data[3]))
			ret.append(utils.sphere(shift+scale*pos,scale*float(data[4]),material=int(data[5]),**kw))

		else:
			raise RuntimeError("Please, specify a correct format output!");
	return ret

def textClumps(fileName,shift=Vector3.Zero,discretization=0,orientation=Quaternion((0,1,0),0.0),scale=1.0,**kw):
	"""Load clumps-members from file, insert them to the simulation.

	:param str filename: file name
	:param str format: the name of output format. Supported `x_y_z_r`(default), `x_y_z_r_clumpId`
	:param [float,float,float] shift: [X,Y,Z] parameter moves the specimen.
	:param float scale: factor scales the given data.
	:param \*\*kw: (unused keyword arguments) is passed to :yref:`sudodem.utils.sphere`
	:returns: list of spheres.

	Lines starting with # are skipped
	"""
	infile = open(fileName,"r")
	lines = infile.readlines()
	infile.close()
	ret=[]

	curClump=[]
	newClumpId = -1

	for line in lines:
		data = line.split()
		if (data[0][0] == "#"): continue
		pos = orientation*Vector3(float(data[0]),float(data[1]),float(data[2]))

		if (newClumpId<0 or newClumpId==int(data[4])):
			idD = curClump.append(utils.sphere(shift+scale*pos,scale*float(data[3]),**kw))
			newClumpId = int(data[4])
		else:
			newClumpId = int(data[4])
			ret.append(O.bodies.appendClumped(curClump,discretization=discretization))
			curClump=[]
			idD = curClump.append(utils.sphere(shift+scale*pos,scale*float(data[3]),**kw))

	if (len(curClump)<>0):
		ret.append(O.bodies.appendClumped(curClump,discretization=discretization))

	# Set the mask to a clump the same as the first member of it
	for i in range(len(ret)):
		O.bodies[ret[i][0]].mask = O.bodies[ret[i][1][0]].mask
	return ret

def text(fileName,shift=Vector3.Zero,scale=1.0,**kw):
	"""Load sphere coordinates from file, returns a list of corresponding bodies; that may be inserted to the simulation with O.bodies.append().

	:param string filename: file which has 4 colums [x, y, z, radius].
	:param [float,float,float] shift: [X,Y,Z] parameter moves the specimen.
	:param float scale: factor scales the given data.
	:param \*\*kw: (unused keyword arguments)	is passed to :yref:`sudodem.utils.sphere`
	:returns: list of spheres.

	Lines starting with # are skipped
	"""

	return textExt(fileName=fileName,format='x_y_z_r',shift=shift,scale=scale,**kw)

def stl(file, dynamic=None,fixed=True,wire=True,color=None,highlight=False,noBound=False,material=-1):
	""" Import geometry from stl file, return list of created facets."""
	imp = STLImporter()
	facets=imp.ymport(file)
	for b in facets:
		b.shape.color=color if color else utils.randomColor()
		b.shape.wire=wire
		b.shape.highlight=highlight
		pos=b.state.pos
		utils._commonBodySetup(b,0,Vector3(0,0,0),material=material,pos=pos,noBound=noBound,dynamic=dynamic,fixed=fixed)
		b.aspherical=False
	return facets

def gts(meshfile,shift=(0,0,0),scale=1.0,**kw):
	""" Read given meshfile in gts format.

	:Parameters:
		`meshfile`: string
			name of the input file.
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`sudodem.utils.facet`
	:Returns: list of facets.
	"""
	import gts,sudodem.pack
	surf=gts.read(open(meshfile))
	surf.scale(scale)
	surf.translate(shift)
	sudodem.pack.gtsSurface2Facets(surf,**kw)

def gmsh(meshfile="file.mesh",shift=Vector3.Zero,scale=1.0,orientation=Quaternion((0,1,0),0.0),**kw):
	""" Imports geometry from mesh file and creates facets.

	:Parameters:
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`orientation`: quaternion
			orientation of the imported mesh
		`**kw`: (unused keyword arguments)
				is passed to :yref:`sudodem.utils.facet`
	:Returns: list of facets forming the specimen.

	mesh files can be easily created with `GMSH <http://www.geuz.org/gmsh/>`_.
	Example added to :ysrc:`examples/regular-sphere-pack/regular-sphere-pack.py`

	Additional examples of mesh-files can be downloaded from
	http://www-roc.inria.fr/gamma/download/download.php
	"""
	infile = open(meshfile,"r")
	lines = infile.readlines()
	infile.close()

	nodelistVector3=[]
	findVerticesString=0

	while (lines[findVerticesString].split()[0]<>'Vertices'): #Find the string with the number of Vertices
		findVerticesString+=1
	findVerticesString+=1
	numNodes = int(lines[findVerticesString].split()[0])

	for i in range(numNodes):
		nodelistVector3.append(Vector3(0.0,0.0,0.0))
	id = 0

	for line in lines[findVerticesString+1:numNodes+findVerticesString+1]:
		data = line.split()
		nodelistVector3[id] = orientation*Vector3(float(data[0])*scale,float(data[1])*scale,float(data[2])*scale)+shift
		id += 1


	findTriangleString=findVerticesString+numNodes
	while (lines[findTriangleString].split()[0]<>'Triangles'): #Find the string with the number of Triangles
		findTriangleString+=1
	findTriangleString+=1
	numTriangles = int(lines[findTriangleString].split()[0])

	triList = []
	for i in range(numTriangles):
		triList.append([0,0,0,0])

	tid = 0
	for line in lines[findTriangleString+1:findTriangleString+numTriangles+1]:
		data = line.split()
		id1 = int(data[0])-1
		id2 = int(data[1])-1
		id3 = int(data[2])-1
		triList[tid][0] = tid
		triList[tid][1] = id1
		triList[tid][2] = id2
		triList[tid][3] = id3
		tid += 1
		ret=[]
	for i in triList:
		a=nodelistVector3[i[1]]
		b=nodelistVector3[i[2]]
		c=nodelistVector3[i[3]]
		ret.append(utils.facet((nodelistVector3[i[1]],nodelistVector3[i[2]],nodelistVector3[i[3]]),**kw))
	return ret

def gengeoFile(fileName="file.geo",shift=Vector3.Zero,scale=1.0,orientation=Quaternion((0,1,0),0.0),**kw):
	""" Imports geometry from LSMGenGeo .geo file and creates spheres.
	Since 2012 the package is available in Debian/Ubuntu and known as python-demgengeo
	http://packages.qa.debian.org/p/python-demgengeo.html

	:Parameters:
		`filename`: string
			file which has 4 colums [x, y, z, radius].
		`shift`: Vector3
			Vector3(X,Y,Z) parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`orientation`: quaternion
			orientation of the imported geometry
		`**kw`: (unused keyword arguments)
				is passed to :yref:`sudodem.utils.sphere`
	:Returns: list of spheres.

	LSMGenGeo library allows one to create pack of spheres
	with given [Rmin:Rmax] with null stress inside the specimen.
	Can be useful for Mining Rock simulation.

	Example: :ysrc:`examples/packs/packs.py`, usage of LSMGenGeo library in :ysrc:`examples/test/genCylLSM.py`.

	* https://answers.launchpad.net/esys-particle/+faq/877
	* http://www.access.edu.au/lsmgengeo_python_doc/current/pythonapi/html/GenGeo-module.html
	* https://svn.esscc.uq.edu.au/svn/esys3/lsm/contrib/LSMGenGeo/"""
	from sudodem.utils import sphere

	infile = open(fileName,"r")
	lines = infile.readlines()
	infile.close()

	numSpheres = int(lines[6].split()[0])
	ret=[]
	for line in lines[7:numSpheres+7]:
		data = line.split()
		pos = orientation*Vector3(float(data[0]),float(data[1]),float(data[2]))
		ret.append(utils.sphere(shift+scale*pos,scale*float(data[3]),**kw))
	return ret

def gengeo(mntable,shift=Vector3.Zero,scale=1.0,**kw):
	""" Imports geometry from LSMGenGeo library and creates spheres.
	Since 2012 the package is available in Debian/Ubuntu and known as python-demgengeo
	http://packages.qa.debian.org/p/python-demgengeo.html

	:Parameters:
		`mntable`: mntable
			object, which creates by LSMGenGeo library, see example
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`sudodem.utils.sphere`

	LSMGenGeo library allows one to create pack of spheres
	with given [Rmin:Rmax] with null stress inside the specimen.
	Can be useful for Mining Rock simulation.

	Example: :ysrc:`examples/packs/packs.py`, usage of LSMGenGeo library in :ysrc:`examples/test/genCylLSM.py`.

	* https://answers.launchpad.net/esys-particle/+faq/877
	* http://www.access.edu.au/lsmgengeo_python_doc/current/pythonapi/html/GenGeo-module.html
	* https://svn.esscc.uq.edu.au/svn/esys3/lsm/contrib/LSMGenGeo/"""
	try:
		from GenGeo import MNTable3D,Sphere
	except ImportError:
		from gengeo import MNTable3D,Sphere
	ret=[]
	sphereList=mntable.getSphereListFromGroup(0)
	for i in range(0, len(sphereList)):
		r=sphereList[i].Radius()
		c=sphereList[i].Centre()
		ret.append(utils.sphere([shift[0]+scale*float(c.X()),shift[1]+scale*float(c.Y()),shift[2]+scale*float(c.Z())],scale*float(r),**kw))
	return ret



def unv(fileName,shift=(0,0,0),scale=1.0,returnConnectivityTable=False,**kw):
	""" Import geometry from unv file, return list of created facets.

		:param string fileName: name of unv file
		:param (float,float,float)|Vector3 shift: (X,Y,Z) parameter moves the specimen.
		:param float scale: factor scales the given data.
		:param \*\*kw: (unused keyword arguments) is passed to :yref:`sudodem.utils.facet`
		:param bool returnConnectivityTable: if True, apart from facets returns also nodes (list of (x,y,z) nodes coordinates) and elements (list of (id1,id2,id3) element nodes ids). If False (default), returns only facets

	unv files are mainly used for FEM analyses (are used by `OOFEM <http://www.oofem.org/>`_ and `Abaqus <http://www.simulia.com/products/abaqus_fea.html>`_), but triangular elements can be imported as facets.
	These files cen be created e.g. with open-source free software `Salome <http://salome-platform.org>`_.

	Example: :ysrc:`examples/test/unv-read/unvRead.py`."""

	class UNVReader:
		# class used in ymport.unv function
		# reads and evaluate given unv file and extracts all triangles
		# can be extended to read tetrahedrons as well
		def __init__(self,fileName,shift=(0,0,0),scale=1.0,returnConnectivityTable=False,**kw):
			self.shift = shift
			self.scale = scale
			self.unvFile = open(fileName,'r')
			self.flag = 0
			self.line = self.unvFile.readline()
			self.lineSplit = self.line.split()
			self.nodes = []
			self.elements = []
			self.read(**kw)
		def readLine(self):
			self.line = self.unvFile.readline()
			self.lineSplit = self.line.split()
		def read(self,**kw):
			while self.line:
				self.evalLine()
				self.line = self.unvFile.readline()
			self.unvFile.close()
			self.createFacets(**kw)
		def evalLine(self):
			self.lineSplit = self.line.split()
			if len(self.lineSplit) <= 1: # eval special unv format
				if   self.lineSplit[0] == '-1':   pass
				elif self.lineSplit[0] == '2411': self.flag = 1; # nodes
				elif self.lineSplit[0] == '2412': self.flag = 2; # edges (lines)
				else: self.flag = 4; # volume elements or other, not interesting for us (at least yet)
			elif self.flag == 1: self.evalNodes()
			elif self.flag == 2: self.evalEdge()
			elif self.flag == 3: self.evalFacet()
			#elif self.flag == 4: self.evalGroup()
		def evalNodes(self):
			self.readLine()
			self.nodes.append((
				self.shift[0]+self.scale*float(self.lineSplit[0]),
				self.shift[1]+self.scale*float(self.lineSplit[1]),
				self.shift[2]+self.scale*float(self.lineSplit[2])))
		def evalEdge(self):
			if self.lineSplit[1]=='41':
				self.flag = 3
				self.evalFacet()
			else:
				self.readLine()
				self.readLine()
		def evalFacet(self):
			if self.lineSplit[1]=='41': # triangle
				self.readLine()
				self.elements.append((
					int(self.lineSplit[0])-1,
					int(self.lineSplit[1])-1,
					int(self.lineSplit[2])-1))
			else: # is not triangle
				self.readLine()
				self.flag = 4
				# can be added function to handle tetrahedrons
		def createFacets(self,**kw):
			self.facets = [utils.facet(tuple(self.nodes[i] for i in e),**kw) for e in self.elements]
	#
	unvReader = UNVReader(fileName,shift,scale,returnConnectivityTable,**kw)
	if returnConnectivityTable:
		return unvReader.facets, unvReader.nodes, unvReader.elements
	return facets


def iges(fileName,shift=(0,0,0),scale=1.0,returnConnectivityTable=False,**kw):
	""" Import triangular mesh from .igs file, return list of created facets.

		:param string fileName: name of iges file
		:param (float,float,float)|Vector3 shift: (X,Y,Z) parameter moves the specimen.
		:param float scale: factor scales the given data.
		:param \*\*kw: (unused keyword arguments) is passed to :yref:`sudodem.utils.facet`
		:param bool returnConnectivityTable: if True, apart from facets returns also nodes (list of (x,y,z) nodes coordinates) and elements (list of (id1,id2,id3) element nodes ids). If False (default), returns only facets
	"""
	nodes,elems = [],[]
	f = open(fileName)
	for line in f:
		if line.startswith('134,'): # read nodes coordinates
			ls = line.split(',')
			v = Vector3(
				float(ls[1])*scale + shift[0],
				float(ls[2])*scale + shift[1],
				float(ls[3])*scale + shift[2]
			)
			nodes.append(v)
		if line.startswith('136,'): # read elements
			ls = line.split(',')
			i1,i2,i3 = int(ls[3])/2, int(ls[4])/2, int(ls[5])/2 # the numbering of nodes is 1,3,5,7,..., hence this int(ls[*])/2
			elems.append( (i1,i2,i3) )
	facets = [utils.facet( ( nodes[e[0]], nodes[e[1]], nodes[e[2]] ), **kw) for e in elems]
	if returnConnectivityTable:
		return facets, nodes, elems
	return facets
