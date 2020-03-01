#pragma once
#ifndef SPHERE_H
#define SPHERE_H
#include<sudodem/core/Shape.hpp>
#include <sudodem/pkg/common/Dispatching.hpp>
#include <sudodem/pkg/common/Aabb.hpp>
#include<sudodem/pkg/common/GLDrawFunctors.hpp>
// HACK to work around https://bugs.launchpad.net/sudodem/+bug/528509
// see comments there for explanation
namespace sudodem{

class Disk: public Shape{
	public:
		Disk(Real _radius): radius(_radius),ref_radius(_radius){createIndex(); }
		virtual ~Disk () {};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(Disk,Shape,"Geometry of spherical particle.",
		((Real,radius,NaN,,"Radius [m]"))
		((Real,ref_radius,NaN,,"reference radius [m]")),
		createIndex(); /*ctor*/
	);
	REGISTER_CLASS_INDEX(Disk,Shape);
};

}
// necessary
using namespace sudodem;

// must be outside sudodem namespace
REGISTER_SERIALIZABLE(Disk);

class Bo1_Disk_Aabb : public BoundFunctor
{
	public :
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se2r&, const Body*);
	FUNCTOR1D(Disk);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Bo1_Disk_Aabb,BoundFunctor,"Functor creating :yref:`Aabb` from :yref:`Disk`.",
		((Real,aabbEnlargeFactor,((void)"deactivated",-1),,"Relative enlargement of the bounding box; deactivated if negative.\n\n.. note::\n\tThis attribute is used to create distant interaction, but is only meaningful with an :yref:`IGeomFunctor` which will not simply discard such interactions: :yref:`Ig2_Disk_Disk_ScGeom::interactionDetectionFactor` should have the same value as :yref:`aabbEnlargeFactor<Bo1_Disk_Aabb::aabbEnlargeFactor>`."))
	);
};

REGISTER_SERIALIZABLE(Bo1_Disk_Aabb);

#ifdef SUDODEM_OPENGL
class Gl1_Disk : public GlShapeFunctor{
	private:
		// for stripes
		static int glStripedDiskList;
		static int glGlutDiskList;
		//Generate GlList for GLUT disk
		void initGlutGlList();
		//Generate GlList for sliced disks
		void initStripedGlList();
		//for regenerating glutDisk list if needed
		static int preSlices;
	public:
		virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
	SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_Disk,GlShapeFunctor,"Renders :yref:`Disk` object",
		((bool,wire,false,,"Only show wireframe (controlled by ``glutSlices`` and ``glutStacks``."))
		((int,Slices,12,(Attr::noSave | Attr::readonly),"Base number of disk slices, multiplied by :yref:`Gl1_Disk::quality` before use); not used with ``stripes`` (see `glut{Solid,Wire}Disk reference <http://www.opengl.org/documentation/specs/glut/spec3/node81.html>`_)"))
	);
	RENDERS(Disk);
};

REGISTER_SERIALIZABLE(Gl1_Disk);
#endif
#endif
