/*************************************************************************
*  Copyright (C) 2017 by Sway Zhao                                       *
*  zhswee@gmail.com                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once

#include<sudodem/pkg/common/GLDrawFunctors.hpp>
#include <sudodem/pkg/common/Dispatching.hpp>
#include<sudodem/core/Scene.hpp>
#include <sudodem/pkg/common/Aabb.hpp>
#include<sudodem/core/Shape.hpp>
#include<sudodem/core/Body.hpp>
#include<sudodem/pkg/fem/Node.hpp>
#include<sudodem/pkg/fem/FEContainer.hpp>

#include<sudodem/pkg/common/Sphere.hpp>

// define this to have topology information about facets enabled;
// it is necessary for TriElementTopologyAnalyzer
// #define FACET_TOPO
typedef Eigen::Matrix<Real,9,1> Vector9r;


// mimick expectation macros that linux has (see e.g. http://kerneltrap.org/node/4705)
#ifndef likely
	#define likely(x) __builtin_expect(!!(x),1)
#endif
#ifndef unlikely
	#define unlikely(x) __builtin_expect(!!(x),0)
#endif//should be moved to Math.hpp


#define MEMBRANE_CONDENSE_DKT

class TriElement : public Shape{
    public:

	virtual ~TriElement();
    //Three nodes
    //Node::id_t nids[3];//node ids in the node container
	int numNodes() const { return nodes.size()==3; }

	/// TriElement's normal
	//Vector3r nf;
	/// Normals of edges
	Vector3r ne[3];
	/// Inscribing cirle radius
	Real icr;
	/// Length of the vertice vectors
	Real vl[3];
	/// Unit vertice vectors
	Vector3r vu[3];
    //gemetry-related func
    Matrix3r triangleInertia(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2);
    Real triangleArea(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2);
    Real TetraVolume(const Vector3r& v);//volume of a tetrhedron constructed with a side point v.
    void lumpMassInertia(Real density);
    //
    Vector3r getCentroid() const;
    Vector3r getNormal() const;
    // closest point on the facet
	Vector3r getNearestPt(const Vector3r& pt);
	vector<Vector3r> outerEdgeNormals() const;
	Real getArea() const;
	Real getPerimeterSq() const;
    // return velocity which is linearly interpolated between velocities of facet nodes, and also angular velocity at that point
	// takes fakeVel in account, with the NaN special case as documented
	std::tuple<Vector3r,Vector3r> interpolatePtLinAngVel(const Vector3r& x) const;
	std::tuple<Vector3r,Vector3r,Vector3r> getOuterVectors() const;
	bool isPointInTriangle(const Vector3r& pt);
	// generic routine: return nearest point on triangle closes to *pt*, given triangle vertices and its normal
	static Vector3r getNearestTrianglePt(const Vector3r& pt, const Vector3r& A, const Vector3r& B, const Vector3r& C, const Vector3r& normal);
    //
    Vector3r closestSegmentPt(const Vector3r& P, const Vector3r& A, const Vector3r& B, Real* normPos);
    Vector3r closestSegmentPt(const Vector3r& P, const Vector3r& A, const Vector3r& B);
    Vector3r glob2loc(const Vector3r& p){ return ori.conjugate()*(p-pos); }
	//
	Vector3r loc2glob(const Vector3r& p){ return ori*p+pos; }
    //
	bool hasRefConf() const { return refRot.size()==3; }
	bool hasBending(){ return KKdkt.size()>0; }
	void pyReset(){ refRot.clear(); }
	void setRefConf(); // use the current configuration as the referential one
	void ensureStiffnessMatrices(const Real& young, const Real& nu, const Real& thickness, bool bending, const Real& bendThickness);
	void updateNode(); // update local coordinate system
	void computeNodalDisplacements(Real dt, bool rotIncr);
	// called by functors to initialize (if necessary) and update
	void stepUpdate(Real dt, bool rotIncr);
	// called from DynDt for updating internal stiffness for given node
	void addIntraStiffnesses(const shared_ptr<Node>&, Vector3r& ktrans, Vector3r& krot) const;

    void applyNodalForces(const Scene* scene,Real dt);
    //
	void postLoad(TriElement&);
    shared_ptr<Node> node1,node2,node3;
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(TriElement,Shape,"TriElement (triangular particle) geometry.",
    ((vector<shared_ptr<Node>>,nodes,,,":yref:`Node` instance associated with this Element."))
    //((shared_ptr<Node>,node1,,,":yref:`Node` instance associated with this Element."))
    //((shared_ptr<Node>,node2,,,":yref:`Node` instance associated with this Element."))
    //((shared_ptr<Node>,node3,,,":yref:`Node` instance associated with this Element."))
    //((vector<Node::id_t>,nids,vector<Node::id_t>(3,Node::ID_NONE),(Attr::triggerPostLoad | Attr::noResize),"Vertex positions in local coordinates."))

    ((Quaternionr,ori,Quaternionr::Identity(),,"Orientation :math:`q` of this node."))
    ((vector<Vector3r>,vertices,vector<Vector3r>(3,Vector3r(NaN,NaN,NaN)),(Attr::triggerPostLoad | Attr::noResize),"Vertex positions in local coordinates."))
    ((Vector3r,normal,Vector3r(NaN,NaN,NaN),(Attr::readonly | Attr::noSave),"TriElement's normal (in local coordinate system)"))
    ((Real,area,NaN,(Attr::readonly | Attr::noSave),"TriElement's area"))((Real,halfThick,0.,,"Geometric thickness (added in all directions)"))
		#ifdef FACET_TOPO
		((vector<Body::id_t>,edgeAdjIds,vector<Body::id_t>(3,Body::ID_NONE),,"TriElement id's that are adjacent to respective edges [experimental]"))
		((vector<Real>,edgeAdjHalfAngle,vector<Real>(3,0),,"half angle between normals of this facet and the adjacent facet [experimental]"))
		#endif
    ((Vector3r,pos,Vector3r::Zero(),,"Position."))
    ((vector<Quaternionr>,refRot,,Attr::readonly,"Rotation applied to nodes to obtain the local coordinate system, computed in the reference configuration. "))
    ((Vector6r,refPos,Vector6r::Zero(),Attr::readonly,"Nodal coordinates in the local coordinate system, in the reference configuration"))
    ((Vector6r,uXy,Vector6r::Zero(),Attr::readonly,"Nodal displacements, stored as ux0, uy0, ux1, uy1, ux1, uy2."))
    ((Real,surfLoad,0.,,"Normal load applied to this facet (positive in the direction of the local normal); this value is multiplied by the current facet's area and equally distributed to nodes."))
    ((Vector6r,phiXy,Vector6r::Zero(),Attr::readonly,"Nodal rotations, only including in-plane rotations (drilling DOF not yet implemented)"))
    ((MatrixXr,KKcst,,,"Stiffness matrix of the element (assembled from the reference configuration when needed for the first time)"))
    ((MatrixXr,KKdkt,,,"Bending stiffness matrix of the element (assembled from the reference configuration when needed for the first time)."))

    /* ((MatrixXr,EBcst,,,"Displacement-stress matrix, for computation of stress tensor in post-processing only."))((MatrixXr,EBdkt,,,"Displacement-stress matrix, for computation of stress tensor in post-processing only.")) */
    #ifdef MEMBRANE_DEBUG_ROT
        ((Vector3r,drill,Vector3r::Zero(),Attr::readonly,"Dirilling rotation (debugging only)"))
			((vector<Quaternionr>,currRot,,Attr::readonly,"What would be the current value of refRot (debugging only!)"))
			((VectorXr,uDkt,,,"DKT displacement vector (saved for debugging only)"))
	#endif
    ,/* ctor */ createIndex();
    ,/*py*/
    .def("setRefConf",&TriElement::setRefConf,"Set the current configuration as the reference one.")
    .def("update",&TriElement::stepUpdate,(py::arg("dt"),py::arg("rotIncr")=false),"Update current configuration; create reference configuration if it does not exist.")
    .def("reset",&TriElement::pyReset,"Reset reference configuration; this forces using the current config as reference when :obj:`update` is called again.")
	);
	DECLARE_LOGGER;
	REGISTER_CLASS_INDEX(TriElement,Shape);
};
REGISTER_SERIALIZABLE(TriElement);
////////
class Bo1_TriElement_Aabb : public BoundFunctor{
	public:
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*);
	FUNCTOR1D(TriElement);
	SUDODEM_CLASS_BASE_DOC(Bo1_TriElement_Aabb,BoundFunctor,"Creates/updates an :yref:`Aabb` of a :yref:`TriElement`.");
};
REGISTER_SERIALIZABLE(Bo1_TriElement_Aabb);
////////
class Gl1_TriElement : public GlShapeFunctor
{
	public:
		virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
	RENDERS(TriElement);
	SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_TriElement,GlShapeFunctor,"Renders :yref:`TriElement` object",
        ((bool,wire,true,,"Only show wireframe"))
		((bool,normals,false,,"In wire mode, render normals of facets and edges; facet's :yref:`colors<Shape::color>` are disregarded in that case."))
	);
};

REGISTER_SERIALIZABLE(Gl1_TriElement);


