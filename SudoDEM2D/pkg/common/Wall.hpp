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
		((int,axis,0,,"Axis of the normal; can be 0,1 for +x, +y respectively (Body's orientation is disregarded for walls)")),
		/*ctor*/createIndex();
	);
	REGISTER_CLASS_INDEX(Wall,Shape);
};
REGISTER_SERIALIZABLE(Wall);

/*! Functor for computing axis-aligned bounding box
    from axis-aligned wall. Has no parameters. */
class Bo1_Wall_Aabb: public BoundFunctor{
	public:
		virtual void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se2r& se2, const Body*);
	FUNCTOR1D(Wall);
	SUDODEM_CLASS_BASE_DOC(Bo1_Wall_Aabb,BoundFunctor,"Creates/updates an :yref:`Aabb` of a :yref:`Wall`");
};
REGISTER_SERIALIZABLE(Bo1_Wall_Aabb);

class Fwall : public Shape {
		public:

	virtual ~Fwall();

	// Postprocessed attributes

	/// Fwall's normal
	//Vector3r nf;
	/// Normals of edge
	//Vector2r ne;
	/// Inscribing cirle radius
	//Real icr;
	/// Length of the vertice vector
	Real vl;
	/// vertice vector
	Vector2r vu;//vertice 1 to vertice 2

	void postLoad(Fwall&);

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Fwall,Shape,"Fwall (triangular particle) geometry.",
		((Vector2r,vertex1,Vector2r(NaN,NaN),(Attr::triggerPostLoad),"Vertex positions in local coordinates."))
		((Vector2r,vertex2,Vector2r(NaN,NaN),(Attr::triggerPostLoad),"Vertex positions in local coordinates."))
		((Vector2r,normal,Vector2r(NaN,NaN),(Attr::readonly | Attr::noSave),"Fwall's normal (in local coordinate system)"))
		#ifdef Fwall_TOPO
		((vector<Body::id_t>,edgeAdjIds,vector<Body::id_t>(2,Body::ID_NONE),,"Fwall id's that are adjacent to respective edges [experimental]"))
		//((vector<Real>,edgeAdjHalfAngle,vector<Real>(3,0),,"half angle between normals of this Fwall and the adjacent Fwall [experimental]"))
		#endif
		,
		/* ctor */ createIndex();,
	);
	DECLARE_LOGGER;
	REGISTER_CLASS_INDEX(Fwall,Shape);
};
REGISTER_SERIALIZABLE(Fwall);

class Bo1_Fwall_Aabb : public BoundFunctor{
	public:
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se2r& se2, const Body*);
	FUNCTOR1D(Fwall);
	SUDODEM_CLASS_BASE_DOC(Bo1_Fwall_Aabb,BoundFunctor,"Creates/updates an :yref:`Aabb` of a :yref:`Fwall`.");
};
REGISTER_SERIALIZABLE(Bo1_Fwall_Aabb);

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

	class Gl1_Fwall : public GlShapeFunctor
	{
		public:
			virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
		RENDERS(Fwall);
		SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_Fwall,GlShapeFunctor,"Renders :yref:`Fwall` object",
		);
	};

	REGISTER_SERIALIZABLE(Gl1_Fwall);

#endif
