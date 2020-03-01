#pragma once
/*
 * =====================================================================================
 *
 *       Filename:  PolySuperellipsoid.h
 *
 *    Description:  PolySuperellipsoid have been defined.
 *
 *        Version:  1.0
 *        Created:  07/22/2015 05:16:35 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  zhswee (zhswee@gmail.com)
 *   Organization: South China University of Technology
 *
 * =====================================================================================
 */

#include<vector>
#include<sudodem/core/Shape.hpp>
#include<sudodem/core/IGeom.hpp>
#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/core/Material.hpp>
#include<sudodem/pkg/common/Aabb.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>
#include<sudodem/pkg/common/Wall.hpp>
#include<sudodem/pkg/common/Facet.hpp>
#include<sudodem/lib/base/openmp-accu.hpp>
#include<sudodem/lib/base/Math.hpp>
//#include <math.h>
//#include <cmath>
//#include <stdlib.h>
//#include <iostream>

//#include <Eigen/Dense>

//typedef Eigen::Quaternion<Real> Quaternionr;
//typedef Eigen::Matrix<Real,3,3> Matrix3r;
typedef Eigen::Matrix<Real, 3, 2> Mat32r;
typedef Eigen::Matrix<Real, 2, 3> Mat23r;
//typedef Eigen::Matrix<Real, 1,2> Vector2r;
typedef Eigen::Matrix<Real, 2,2> Matrix2d;//move to Math.hpp later

//using namespace Eigen;
//**********************************************************************************
class PolySuperellipsoid: public Shape{
	public:
		//constructor
		PolySuperellipsoid(Vector2r eps, Vector6r halflen);//{/*createIndex();*/ rxyz = Vector3r(x,y,z); eps = Vector2r(ep1,ep2); Initial();};
		/*PolySuperellipsoid(Vector3r xyz, Vector2r _eps):rxyz(xyz),eps(_eps)
		{
		createIndex();

	        //initialize to calculate some basic geometric parameters.
	        Initial();
	        };*/
		virtual ~PolySuperellipsoid();
		void Initial();				//calculate some basic geometric parameters.
		bool IsInitialized(){return init;}
		/////////get
		//double getrx(){return rx;}
		//double getry(){return ry;}
		//double getrz(){return rz;}
		Vector6r getrxyz(){return rxyz;}
		Vector6r getrxyz_ref(){return rxyz_ref;}
		//double geteps1(){return eps1;}
		//double geteps2(){return eps2;}
		Vector2r geteps(){return eps;}
		double getr_max(){return r_max;}
		double getr_max_ref(){return r_max_ref;}
		Vector3r getPosition(){return Position;}
		double getVolume(){Initial();return Volume;}
		Vector3r getMassCenter(){Initial();return massCenter;}
		Vector3r getMassCenter_ref(){return massCenter_ref;}
		Vector3r getInertia(){Initial();return Inertia;}
		Quaternionr getOrientation(){return Orientation;}
    bool isInside(Vector3r p);

		void setPosition(Vector3r p){Position = p;}
		void setOrientation(Quaternionr Ori){
		Ori.normalize();Orientation = Ori;
		//rotation matrices
    rot_mat2local = (Orientation).conjugate().toRotationMatrix();//to particle's system
	  rot_mat2global = (Orientation).toRotationMatrix();//to global system
		}
		void setRxyz(Vector6r rxyz_in){rxyz = rxyz_in;}
		void setMassCenter(Vector3r mc){massCenter = mc;}
		void setRmax(double rm){r_max = rm;}

    //rotation matrices
    Matrix3r rot_mat2local; //(Orientation).conjugate().toRotationMatrix();//to particle's system
	  Matrix3r rot_mat2global; // (Orientation).toRotationMatrix();//to global system
		Vector3r getSurface(Vector2r phi) const;//at the local with respect to the geometric center
		Vector3r getSurfaceMC(Vector2r phi) const;//at the local with respect to the mass center
		Vector3r getNormal(Vector2r);
		Vector2r Normal2Phi(Vector3r);
		Vector3r Normal2SurfaceMC(Vector3r n) const;
		Vector3r Normal2SurfaceMCgl(Vector3r gn) const;
		Vector3r support(Vector3r n) const;
	protected:
		//double rx, ry, rz, eps1, eps2;//the main parameters of a superquadric
		//r_max:the maximum r in the three axes
	  double r_max, r_max_ref;
	  //sign of performed initialization
		bool init;
		Vector3r Position;			//the center position of a superquadric, i.e (x, y, z)
		double Volume;
		Vector6r rxyz_ref;
		Vector3r massCenter;
		Vector3r massCenter_ref;//used in the expansion method
		Vector3r Inertia;			//the pricipal moments of inetia
		Quaternionr Orientation;	//orientation
		Vector3r Surface;			//surface equations represented by the sureface parameters Phi1 and Phi2, where Phi1 belongs [-pi,pi], Phi2 in [-pi/2., pi/2.].
		Vector3r CurrentP;			//the local current point (x,y,z) that is used to calculate some parameters.
		Vector2r CurrentPhi;		//the local current phi1,phi2
		Vector3r CurrentNorm;		//the local current n1,n2,n3
		Vector3r ContactNorm;		//the vector of contact direction
public:
    int m_GL_slices;//
    int m_glSolidListNum;//glList number
    int m_glWireListNum;//glList number
		SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(PolySuperellipsoid,Shape,"PolySuperellipsoid shape.",
                        //((Vector3r,rxyz,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
                        //((double,rx,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
                        //((double,ry,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
                        //((double,rz,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
												((Vector6r,rxyz,,,"length of half-axis in the six driections, x,-x,y,-y,z,-z, local coordinate system"))
                        ((double,eps1,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
                        ((double,eps2,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
                        ((bool,isSphere,false,,"the shape is non-spherical by default. Using this flag for accelerating spherical systems"))
                        ((Vector2r,eps,,,"length of half-axis in the three driections, x,y,z, local coordinate system")),
			/*init*/,
			/*ctor*/
			createIndex();
			init = false;
            m_glSolidListNum = -1;
            m_glWireListNum = -1;
			m_GL_slices = 10,
			.def("Initial",&PolySuperellipsoid::Initial,"Initialization")
			.def("getVolume",&PolySuperellipsoid::getVolume,"return particle's volume")
			.def("getMassCenter",&PolySuperellipsoid::getMassCenter,"mass center at the local coordinate system with a fixed origin at the geometry center")
			.def("getrxyz",&PolySuperellipsoid::getrxyz,"return particle's rxyz")
			//.def("getry",&PolySuperellipsoid::getry,"return particle's ry")
			//.def("getrz",&PolySuperellipsoid::getrz,"return particle's rz")
			.def("geteps",&PolySuperellipsoid::geteps,"return particle's eps1 and eps2")
			//.def("geteps2",&PolySuperellipsoid::geteps2,"return particle's eps2")
			.def("getInertia",&PolySuperellipsoid::getInertia,"return particle's inertia tensor")
			.def("getOri",&PolySuperellipsoid::getOrientation,"return particle's orientation")
			.def("getCentroid",&PolySuperellipsoid::getPosition,"return particle's centroid")
		);
		REGISTER_CLASS_INDEX(PolySuperellipsoid,Shape);
};
REGISTER_SERIALIZABLE(PolySuperellipsoid);


//***************************************************************************
/*! Collision configuration for PolySuperellipsoid and something.
 * This is expressed as penetration volume properties: centroid, volume, depth ...
 *
 * Self-contained. */
class PolySuperellipsoidGeom: public IGeom{
	public:
		virtual ~PolySuperellipsoidGeom();

		//precompute data for shear evaluation
		void precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector3r& currentNormal, bool isNew, const Vector3r& shift2);
		Vector3r& rotate(Vector3r& shearForce) const;
		//sep_plane is a code storing plane, that previously separated two polyhedras. It is used for faster detection of non-overlap.
		//std::vector<int> sep_plane;
		//bool isShearNew;
	//protected:
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(PolySuperellipsoidGeom,IGeom,"Geometry of interaction between 2 :yref:`vector<PolySuperellipsoid>`, including volumetric characteristics",

			((Real,PenetrationDepth,0.,,"PenetrationDepth"))
			((Vector2r,contactAngle,Vector2r::Zero(),,"Normal direction (in local coordinates)"))
			((Vector3r,contactPoint,Vector3r::Zero(),,"Contact point (global coords), center of the PenetrationDepth"))
			((Vector3r,shearInc,Vector3r::Zero(),,"Shear displacement increment in the last step"))
			((Vector3r,relativeVn,Vector3r::Zero(),,"relative velocity in the normal at the contact"))
			((Vector3r,relativeVs,Vector3r::Zero(),,"relative velocity in the tangential at the contact"))
			((Vector3r,normal,Vector3r::Zero(),,"Contact normal"))
			((Vector3r,twist_axis,Vector3r::Zero(),,""))
			((Vector3r,point1,Vector3r::Zero(),,"point 1 of the closest two points."))
			((Vector3r,point2,Vector3r::Zero(),,"point 2 of the closest two points."))
			((Vector3r,point11,Vector3r::Zero(),,"direction from point 1.FOR TESTING"))
			((Vector3r,point22,Vector3r::Zero(),,"direction from point 2.FOR TESTING"))
			((bool,isShearNew,true,,""))
			((Vector3r,orthonormal_axis,Vector3r::Zero(),,"")),
			createIndex();
			//sep_plane.assign(3,0);
		);
		//FUNCTOR2D(Tetra,Tetra);
		REGISTER_CLASS_INDEX(PolySuperellipsoidGeom,IGeom);
};
REGISTER_SERIALIZABLE(PolySuperellipsoidGeom);
//***************************************************************************
/*! Creates Aabb from PolySuperellipsoid.
 *
 * Self-contained. */
class Bo1_PolySuperellipsoid_Aabb: public BoundFunctor{
	public:
		void go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se3r& se3, const Body*);
		FUNCTOR1D(PolySuperellipsoid);
		SUDODEM_CLASS_BASE_DOC(Bo1_PolySuperellipsoid_Aabb,BoundFunctor,"Create/update :yref:`Aabb` of a :yref:`PolySuperellipsoid`");
};
REGISTER_SERIALIZABLE(Bo1_PolySuperellipsoid_Aabb);

//***************************************************************************
/*! Elastic material */
class PolySuperellipsoidMat: public Material{
	public:
		 PolySuperellipsoidMat(double N, double S, double F){Kn=N; Ks=S; frictionAngle=F;};
		 double GetStrength(){return strength;};
	virtual ~PolySuperellipsoidMat(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(PolySuperellipsoidMat,Material,"Elastic material with Coulomb friction.",
		((Real,Kn,1e8,,"Normal stiffness (N/m)."))
		((Real,Ks,1e5,,"Shear stiffness (N/m)."))
		((Real,frictionAngle,.5,,"Contact friction angle (in radians)."))
		((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
		((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
		((bool,IsSplitable,0,,"To be splitted ... or not"))
		((double,strength,100,,"Stress at whis polyhedra of volume 4/3*pi [mm] breaks.")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(PolySuperellipsoidMat,Material);
};
REGISTER_SERIALIZABLE(PolySuperellipsoidMat);
//***************************************************************************
//use FrictPhys as the basic class
class PolySuperellipsoidPhys: public FrictPhys{
	public:
	virtual ~PolySuperellipsoidPhys(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(PolySuperellipsoidPhys,FrictPhys,"Simple elastic material with friction for volumetric constitutive laws",
		//((Real,kn,0,,"Normal stiffness"))
		//((Vector3r,normalForce,Vector3r::Zero(),,"Normal force after previous step (in global coordinates)."))
		//((Real,ks,0,,"Shear stiffness"))
		//((Vector3r,shearForce,Vector3r::Zero(),,"Shear force after previous step (in global coordinates)."))
		// Contact damping ratio as for linear elastic contact law
		((Vector3r,normalViscous,Vector3r::Zero(),,"Normal viscous component"))
		((Vector3r,shearViscous,Vector3r::Zero(),,"Shear viscous component"))
		((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
		((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$.")),
		//((Real,tangensOfFrictionAngle,0.,,"tangens of angle of internal friction")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(PolySuperellipsoidPhys,FrictPhys);
};
REGISTER_SERIALIZABLE(PolySuperellipsoidPhys);
//***************************************************************************
#ifdef SUDODEM_OPENGL
	#include<sudodem/pkg/common/GLDrawFunctors.hpp>
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	#include<sudodem/lib/opengl/GLUtils.hpp>
	#include<GL/glu.h>
	#include<sudodem/pkg/dem/Shop.hpp>

	/*! Draw PolySuperellipsoid using OpenGL */
	class Gl1_PolySuperellipsoid: public GlShapeFunctor{
		private:
			//static vector<Vector3r> vertices;
            static int glSolidList;
		    static int glWireList;
            void initSolidGlList(PolySuperellipsoid*);
            void initWireGlList(PolySuperellipsoid*);
            static int pre_slices;
			void drawSlice(Vector3r &p1,Vector3r &p2,Vector3r &p3,Vector3r &p4);

		public:
			virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
			SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_PolySuperellipsoid,GlShapeFunctor,"Renders :yref:`PolySuperellipsoid` object",
			((bool,wire,true,,"Only show wireframe"))
			((int,slices,5,,"Base number of slices."))
			);
			RENDERS(PolySuperellipsoid);
	};
	REGISTER_SERIALIZABLE(Gl1_PolySuperellipsoid);

	struct Gl1_PolySuperellipsoidGeom: public GlIGeomFunctor{
		RENDERS(PolySuperellipsoidGeom);
		void go(const shared_ptr<IGeom>&, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool);
		void draw(const shared_ptr<IGeom>&);
		SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_PolySuperellipsoidGeom,GlIGeomFunctor,"Render :yref:`PolySuperellipsoidGeom` geometry.",
		);
	};
	REGISTER_SERIALIZABLE(Gl1_PolySuperellipsoidGeom);
#endif


//***************************************************************************
class Ip2_PolySuperellipsoidMat_PolySuperellipsoidMat_PolySuperellipsoidPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
			const shared_ptr<Material>& b2,
			const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(PolySuperellipsoidMat,PolySuperellipsoidMat);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_PolySuperellipsoidMat_PolySuperellipsoidMat_PolySuperellipsoidPhys,IPhysFunctor,"",
	//((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))//will crash:Abort (core dumped)
	//((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
	);
};
REGISTER_SERIALIZABLE(Ip2_PolySuperellipsoidMat_PolySuperellipsoidMat_PolySuperellipsoidPhys);
//***************************************************************************
/*! Calculate physical response based on penetration configuration given by TTetraGeom. */

class PolySuperellipsoidLaw: public LawFunctor{
	OpenMPAccumulator<Real> plasticDissipation;
	virtual bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
	Real elasticEnergy ();
	Real getPlasticDissipation();
	void initPlasticDissipation(Real initVal=0);
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(PolySuperellipsoidLaw,LawFunctor,"Calculate physical response of 2 :yref:`vector<PolySuperellipsoid>` in interaction, based on penetration configuration given by :yref:`PolySuperellipsoidGeom`.",
	((Vector3r,shearForce,Vector3r::Zero(),,"Shear force from last step"))
	((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts. This will trace only plastic energy in this law, see O.trackEnergy for a more complete energies tracing"))
	((int,plastDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for plastic dissipation (with O.trackEnergy)"))
	((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
	,,
	.def("elasticEnergy",&PolySuperellipsoidLaw::elasticEnergy,"Compute and return the total elastic energy in all \"FrictPhys\" contacts")
	.def("plasticDissipation",&PolySuperellipsoidLaw::getPlasticDissipation,"Total energy dissipated in plastic slips at all FrictPhys contacts. Computed only if :yref:`Law2_ScGeom_FrictPhys_CundallStrack::traceEnergy` is true.")
	.def("initPlasticDissipation",&PolySuperellipsoidLaw::initPlasticDissipation,"Initialize cummulated plastic dissipation to a value (0 by default).")
	);
	FUNCTOR2D(PolySuperellipsoidGeom,PolySuperellipsoidPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PolySuperellipsoidLaw);
//********************************************************************************
//shared_ptr<Body> NewPolySuperellipsoid(double x, double y, double z, double ep1, double ep2, shared_ptr<Material> mat,bool rotate);
