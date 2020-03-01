# coding: utf-8
# 2009 © Václav Šmilauer <eudoxos@arcig.cz>

"""
Functions for accessing sudodem's internals; only used internally.
"""
import sys
from sudodem import wrapper
from sudodem._customConverters import *
from sudodem import runtime
O=wrapper.Omega()

def childClasses(base,recurse=True,includeBase=False):
	"""Enumerate classes deriving from given base (as string), recursively by default. Returns set."""
	ret=set(O.childClassesNonrecursive(base)); ret2=set();
	if includeBase: ret|=set([base])
	if not recurse: return ret
	for bb in ret:
		ret2|=childClasses(bb)
	return ret | ret2

_allSerializables=childClasses('Serializable')
## set of classes for which the proxies were created
_proxiedClasses=set()

## deprecated classes
# if old class name is used, the new object is constructed and a warning is issued about old name being used
# keep chronologically ordered, oldest first; script/rename-class.py appends at the end
_deprecated={
	'GranularMat':'FrictMat', # Sun Jan 10 09:26:45 2010, vaclav@flux
	'SimpleElasticRelationships':'Ip2_FrictMat_FrictMat_NormShearPhys', # Sun Jan 10 09:28:17 2010, vaclav@flux
	'NormalInteraction':'NormPhys', # Sun Jan 10 09:28:56 2010, vaclav@flux
	'NormalShearInteraction':'NormShearPhys', # Sun Jan 10 09:29:22 2010, vaclav@flux
	'ElasticMat':'ElastMat', # Sun Jan 10 09:53:15 2010, vaclav@flux
	'ElasticContactInteraction':'FrictPhys', # Sun Jan 10 09:57:59 2010, vaclav@flux
	'ef2_Spheres_Elastic_ElasticLaw':'Law2_ScGeom_FrictPhys_Basic', # Sun Jan 10 09:59:42 2010, vaclav@flux
	'Ip2_FrictMat_FrictMat_NormShearPhys':'Ip2_FrictMat_FrictMat_FrictPhys', # Sun Jan 10 10:07:40 2010, vaclav@flux
	'GLDrawCpmPhys':'Gl1_CpmPhys', # Sat Feb  6 14:46:08 2010, vaclav@flux
	'RockJointLawRelationships':'Ip2_2xCohFrictMat_NormalInelasticityPhys', # Mon Feb  8 11:17:59 2010, jduriez@c1solimara-l
	'TetraBang':'TTetraGeom', # Tue Feb  9 10:21:24 2010, vaclav@flux
	'TetraMold':'Tetra', # Tue Feb  9 10:22:15 2010, vaclav@flux
	'TetraAABB':'Bo1_Tetra_Aabb', # Tue Feb  9 10:22:33 2010, vaclav@flux
	'Tetra2TetraBang':'Ig2_Tetra_Tetra_TTetraGeom', # Tue Feb  9 10:23:19 2010, vaclav@flux
	'TetraLaw':'TetraVolumetricLaw', # Tue Feb  9 10:24:10 2010, vaclav@flux
	'DirecResearchEngine':'Disp2DPropLoadEngine', # Wed Mar 10 12:23:42 2010, jduriez@c1solimara-l
	'Ip2_BMP_BMP_CSPhys':'Ip2_2xFrictMat_CSPhys', # Wed Mar 10 15:08:56 2010, eudoxos@frigo
	'NormalInelasticityLaw':'Law2_ScGeom_NormalInelasticityPhys_NormalInelasticity', # Wed Mar 17 15:50:59 2010, jduriez@c1solimara-l
	'CapillaryCohesiveLaw':'CapillaryLaw', # Tue Mar 30 14:11:36 2010, sch50p@fluent-ph
	'SimpleElasticRelationshipsWater':'Ip2_Frictmat_FrictMat_CapillaryLawPhys', # Tue Mar 30 14:20:36 2010, sch50p@fluent-ph
	'CapillaryLaw':'Law2_ScGeom_CapillaryPhys_Capillarity', # Wed Mar 31 09:23:36 2010, sch50p@fluent-ph
	'CapillaryParameters':'CapillaryPhys', # Wed Mar 31 09:25:03 2010, sch50p@fluent-ph
	'Ip2_FrictMat_FrictMat_CapillaryLawPhys':'Ip2_FrictMat_FrictMat_CapillaryPhys', # Wed Mar 31 09:26:04 2010, sch50p@fluent-ph
	'SimpleViscoelasticMat':'ViscElMat', # Fri Apr  9 19:25:38 2010, vaclav@flux
	'SimpleViscoelasticPhys':'ViscElPhys', # Fri Apr  9 19:26:34 2010, vaclav@flux
	'Law2_Spheres_Viscoelastic_SimpleViscoelastic':'Law2_ScGeom_ViscElPhys_Basic', # Fri Apr  9 19:28:02 2010, vaclav@flux
	'Ip2_SimleViscoelasticMat_SimpleViscoelasticMat_SimpleViscoelasticPhys':'Ip2_ViscElMat_ViscElMat_ViscElPhys', # Fri Apr  9 19:28:48 2010, vaclav@flux
	'MomentEngine':'TorqueEngine', # Sun May  2 16:09:34 2010, vaclav@flux
	'JumpChangeSe3':'StepDisplacer', # Sun May  2 16:14:21 2010, vaclav@flux
	'Ip2_2xCohFrictMat_NormalInelasticityPhys':'Ip2_2xNormalInelasticMat_NormalInelasticityPhys', # Fri Jun  4 15:36:41 2010, jduriez@c1solimara-l
	'Ip2_2xCohFrictMat_NormalInelasticityPhys':'Ip2_2xNormalInelasticMat_NormalInelasticityPhys', # Fri Jun  4 15:37:01 2010, jduriez@c1solimara-l
	'Ip2_2xCohFrictMat_NormalInelasticityPhys':'Ip2_2xNormalInelasticMat_NormalInelasticityPhys', # Fri Jun  4 15:37:16 2010, jduriez@c1solimara-l
	'OpenGLRenderingEngine':'OpenGLRenderer', # Sat Jul 24 06:04:13 2010, vaclav@flux
	'PeriodicPythonRunner':'PyRunner', # Wed Sep  1 16:41:50 2010, chia@engs-018373
	'InteractionDispatchers':'InteractionLoop', # Mon Sep 27 13:44:54 2010, chia@engs-018373
	'InteractionGeometry':'IGeom', # Thu Sep 30 10:39:24 2010, chia@engs-018373
	'InteractionPhysics':'IPhys', # Thu Sep 30 10:39:43 2010, chia@engs-018373
	'InteractionGeometryFunctor':'IGeomFunctor', # Thu Sep 30 14:27:43 2010, chia@engs-018373
	'InteractionPhysicsFunctor':'IPhysFunctor', # Thu Sep 30 14:27:53 2010, chia@engs-018373
	'InteractionPhysicsDispatcher':'IPhysDispatcher', # Thu Sep 30 14:28:12 2010, chia@engs-018373
	'InteractionGeometryDispatcher':'IGeomDispatcher', # Thu Sep 30 14:28:22 2010, chia@engs-018373
	'Law2_ScGeom_FrictPhys_Basic':'Law2_ScGeom_FrictPhys_CundallStrack', # Wed Oct 13 17:40:42 2010, bchareyre@dt-rv020
	'Law2_ScGeom_CohFrictPhys_ElasticPlastic':'Law2_ScGeom_CohFrictPhys_CohesionMoment', # Wed Oct 13 17:47:09 2010, bchareyre@dt-rv020
	'Law2_ScGeom_NormalInelasticityPhys_NormalInelasticity':'Law2_ScGeom6D_NormalInelasticityPhys_NormalInelasticity', # Thu Oct 28 18:09:16 2010, jduriez@c1solimara-l
	'SpiralEngine':'HelixEngine', # Sat Oct 30 05:52:06 2010, vaclav@flux
	'InterpolatingSpiralEngine':'InterpolatingHelixEngine', # Sat Oct 30 05:52:21 2010, vaclav@flux
	'SpiralInteractionLocator2d':'HelixInteractionLocator2d', # Sat Oct 30 05:53:04 2010, vaclav@flux
	'Law2_ScGeom_CohFrictPhys_CohesionMoment':'Law2_ScGeom6D_CohFrictPhys_CohesionMoment', # Fri Nov 12 18:45:23 2010, bchareyre@dt-rv020
	'Ig2_ChainedCylinder_ChainedCylinder_ScGeom':'Ig2_ChainedCylinder_ChainedCylinder_ScGeom6D', # Fri Nov 12 18:47:44 2010, bchareyre@dt-rv020
	'Ip2_2xCohFrictMat_CohFrictPhys':'Ip2_CohFrictMat_CohFrictMat_CohFrictPhys', # Tue Dec 21 16:11:28 2010, bchareyre@dt-rv020
	'Ig2_Sphere_Sphere_L3Geom_Inc':'Ig2_Sphere_Sphere_L3Geom', # Sun Dec 26 11:11:05 2010, vaclav@flux
	'Ig2_Wall_Sphere_L3Geom_Inc':'Ig2_Wall_Sphere_L3Geom', # Sun Dec 26 11:11:30 2010, vaclav@flux
	'Ig2_Sphere_Sphere_L6Geom_Inc':'Ig2_Sphere_Sphere_L6Geom', # Sun Dec 26 11:11:53 2010, vaclav@flux
	### END_RENAMED_CLASSES_LIST ### (do not delete this line; scripts/rename-class.py uses it
}

def updateScripts(scripts):
	## Thanks goes to http://code.activestate.com/recipes/81330-single-pass-multiple-replace/
	from UserDict import UserDict
	import re,os
	class Xlator(UserDict):
		"An all-in-one multiple string substitution class; adapted to match only whole words"
		def _make_regex(self):
			"Build a regular expression object based on the keys of the current dictionary"
			return re.compile(r"(\b%s\b)" % "|".join(self.keys()))  ## adapted here
		def __call__(self, mo):
			"This handler will be invoked for each regex match"
			# Count substitutions
			self.count += 1 # Look-up string
			return self[mo.string[mo.start():mo.end()]]
		def xlat(self, text):
			"Translate text, returns the modified text."
			# Reset substitution counter
			self.count = 0
			# Process text
			return self._make_regex().sub(self, text)
	# use the _deprecated dictionary for translation, but only when matching on words boundary
	xlator=Xlator(_deprecated)
	if len(scripts)==0: print "No scripts given to --update. Nothing to do."
	for s in scripts:
		if not s.endswith('.py'): raise RuntimeError("Refusing to do --update on file '"+s+"' (not *.py)")
		txt=open(s).read()
		txt2=xlator.xlat(txt)
		if xlator.count==0: print "%s: already up-to-date."%s
		else:
			os.rename(s,s+'~')
			out=open(s,'w'); out.write(txt2); out.close()
			print "%s: %d subtitution%s made, backup in %s~"%(s,xlator.count,'s' if xlator.count>1 else '',s)


def cxxCtorsDict(proxyNamespace=__builtins__):
	"""Return dictionary of class constructors for sudodem's c++ types, which should be used to update a namespace.

	Root classes are those that are directly wrapped by boost::python. These are only put to the dict.

	Derived classes (from these root classes) are faked by creating a callable which invokes appropriate root class constructor with the derived class parameter and passes remaining arguments to it.

	Classes that are neither root nor derived are exposed via callable object that constructs a Serializable of given type and passes the parameters.
	"""
	proxyNamespace={}

	import sudodem.wrapper
	for c in _allSerializables:
		try:
			proxyNamespace[c]=sudodem.wrapper.__dict__[c]
		except KeyError: pass # not registered properly

	# deprecated names
	for oldName in _deprecated.keys():
		class warnWrap:
			def __init__(self,_old,_new):
				# assert(proxyNamespace.has_key(_new))
				self.old,self.new=_old,_new
			def __call__(self,*args,**kw):
				import warnings; warnings.warn("Class `%s' was renamed to (or replaced by) `%s', update your code! (you can run 'sudodem --update script.py' to do that automatically)"%(self.old,self.new),DeprecationWarning,stacklevel=2);
				return sudodem.wrapper.__dict__[self.new](*args,**kw)
		proxyNamespace[oldName]=warnWrap(oldName,_deprecated[oldName])
	return proxyNamespace


# consistency check
# if there are no serializables, then plugins were not loaded yet, probably
if(len(_allSerializables)==0):
	raise ImportError("No classes deriving from Serializable found; you must call sudodem.boot.initialize to load plugins before importing sudodem.system.")
