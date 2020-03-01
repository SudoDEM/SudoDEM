// © 2009 Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/core/Shape.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>


/*! Object representing infinite plane aligned with the coordinate system (axis-aligned wall). */
class Wall: public Shape{
	public:
		virtual ~Wall(); // vtable
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(Wall,Shape,"Object representing infinite plane aligned with the coordinate system (axis-aligned wall).",
		((int,sense,0,,"Which side of the wall interacts: -1 for negative only, 0 for both, +1 for positive only"))
		((int,axis,0,,"Axis of the normal; can be 0,1,2 for +x, +y, +z respectively (Body's orientation is disregarded for walls)")),
		/*ctor*/createIndex();
	);
	REGISTER_CLASS_INDEX(Wall,Shape);
};
REGISTER_SERIALIZABLE(Wall);

/*! Functor for computing axis-aligned bounding box
    from axis-aligned wall. Has no parameters. */
class Bo1_Wall_Aabb: public BoundFunctor{
	public:
		virtual void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*);
	FUNCTOR1D(Wall);
	SUDODEM_CLASS_BASE_DOC(Bo1_Wall_Aabb,BoundFunctor,"Creates/updates an :yref:`Aabb` of a :yref:`Wall`");
};
REGISTER_SERIALIZABLE(Bo1_Wall_Aabb);
#ifdef SUDODEM_OPENGL
	#include<sudodem/pkg/common/GLDrawFunctors.hpp>
	class Gl1_Wall: public GlShapeFunctor{
		public:
			virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
		RENDERS(Wall);
		SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_Wall,GlShapeFunctor,"Renders :yref:`Wall` object",
			((int,div,20,,"Number of divisions of the wall inside visible scene part."))
		);
	};
	REGISTER_SERIALIZABLE(Gl1_Wall);
#endif

