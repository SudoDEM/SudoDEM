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

class Sphere: public Shape{
	public:
		Sphere(Real _radius): radius(_radius){ createIndex(); }
		virtual ~Sphere () {};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(Sphere,Shape,"Geometry of spherical particle.",
		((Real,radius,NaN,,"Radius [m]")),
		createIndex(); /*ctor*/
	);
	REGISTER_CLASS_INDEX(Sphere,Shape);
};

}
// necessary
using namespace sudodem;

// must be outside sudodem namespace
REGISTER_SERIALIZABLE(Sphere);

class Bo1_Sphere_Aabb : public BoundFunctor
{
	public :
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r&, const Body*);
	FUNCTOR1D(Sphere);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Bo1_Sphere_Aabb,BoundFunctor,"Functor creating :yref:`Aabb` from :yref:`Sphere`.",
		((Real,aabbEnlargeFactor,((void)"deactivated",-1),,"Relative enlargement of the bounding box; deactivated if negative.\n\n.. note::\n\tThis attribute is used to create distant interaction, but is only meaningful with an :yref:`IGeomFunctor` which will not simply discard such interactions: :yref:`Ig2_Sphere_Sphere_ScGeom::interactionDetectionFactor` should have the same value as :yref:`aabbEnlargeFactor<Bo1_Sphere_Aabb::aabbEnlargeFactor>`."))
	);
};

REGISTER_SERIALIZABLE(Bo1_Sphere_Aabb);

class Gl1_Sphere : public GlShapeFunctor{
	private:
		// for stripes
		static vector<Vector3r> vertices, faces;
		static int glStripedSphereList;
		static int glGlutSphereList;
		void subdivideTriangle(Vector3r& v1,Vector3r& v2,Vector3r& v3, int depth);
		//Generate GlList for GLUT sphere
		void initGlutGlList();
		//Generate GlList for sliced spheres
		void initStripedGlList();
		//for regenerating glutSphere list if needed
		static Real prevQuality;
	public:
		virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
	SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_Sphere,GlShapeFunctor,"Renders :yref:`Sphere` object",
		((Real,quality,1.0,,"Change discretization level of spheres. quality>1  for better image quality, at the price of more cpu/gpu usage, 0<quality<1 for faster rendering. If mono-color spheres are displayed (:yref:`Gl1_Sphere::stripes` = False), quality mutiplies :yref:`Gl1_Sphere::glutSlices` and :yref:`Gl1_Sphere::glutStacks`. If striped spheres are displayed (:yref:`Gl1_Sphere::stripes` = True), only integer increments are meaningfull : quality=1 and quality=1.9 will give the same result, quality=2 will give finer result."))
		((bool,wire,false,,"Only show wireframe (controlled by ``glutSlices`` and ``glutStacks``."))
		((bool,stripes,false,,"In non-wire rendering, show stripes clearly showing particle rotation."))
		((bool,localSpecView,true,,"Compute specular light in local eye coordinate system."))
		((int,glutSlices,12,(Attr::noSave | Attr::readonly),"Base number of sphere slices, multiplied by :yref:`Gl1_Sphere::quality` before use); not used with ``stripes`` (see `glut{Solid,Wire}Sphere reference <http://www.opengl.org/documentation/specs/glut/spec3/node81.html>`_)"))
		((int,glutStacks,6,(Attr::noSave | Attr::readonly),"Base number of sphere stacks, multiplied by :yref:`Gl1_Sphere::quality` before use; not used with ``stripes`` (see `glut{Solid,Wire}Sphere reference <http://www.opengl.org/documentation/specs/glut/spec3/node81.html>`_)"))
	);
	RENDERS(Sphere);
};

REGISTER_SERIALIZABLE(Gl1_Sphere);
#endif
