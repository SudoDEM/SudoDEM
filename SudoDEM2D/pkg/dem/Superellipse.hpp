#pragma once
/*
 * =====================================================================================
 *
 *       Filename:  superellipse.h
 *
 *    Description:  superellipse have been defined.
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
//#include<sudodem/pkg/common/Facet.hpp>
#include<sudodem/lib/base/openmp-accu.hpp>

//#include <math.h>
//#include <cmath>
//#include <stdlib.h>
//#include <iostream>

//#include <Eigen/Dense>

typedef Eigen::Quaternion<Real> Quaternionr;
typedef Eigen::Matrix<Real,3,3> Matrix3r;
typedef Eigen::Matrix<Real, 3, 2> Mat32r;
typedef Eigen::Matrix<Real, 2, 3> Mat23r;
//typedef Eigen::Matrix<Real, 1,2> Vector2r;
typedef Eigen::Matrix<Real, 2,2> Matrix2d;
//using namespace Eigen;
//**********************************************************************************
class Superellipse: public Shape{
	public:
		//constructor
		Superellipse(double x, double y, double ep);//{/*createIndex();*/ rxy = Vector2r(x,y); eps = Vector2r(ep1,ep2); Initial();};
		/*Superellipse(Vector2r xyz, Vector2r _eps):rxyz(xyz),eps(_eps)
		{
		createIndex();

	        //initialize to calculate some basic geometric parameters.
	        Initial();
	        };*/
		virtual ~Superellipse();
		void Initial();				//calculate some basic geometric parameters.
		bool IsInitialized(){return init;}
		/////////get
		//double getrx(){return rx;}
		//double getry(){return ry;}
		//double getrz(){return rz;}
		Vector2r getrxy(){return rxy;}
		//double geteps1(){return eps1;}
		//double geteps2(){return eps2;}
		Real geteps(){return eps;}
		double getr_max(){return r_max;}
		Vector2r getPosition(){return Position;}
		double getArea(){Initial();return Area;}
		Real getInertia(){Initial();return Inertia;}
		Rotationr getOrientation(){return Orientation;}
    bool isInside(Vector2r p);

		//void setrx(double x){rx = x;}
		//void setry(double y){ry = y;}
		//void setrz(double z){rz = z;}
		//void seteps1(double ep1){eps1 = ep1;}
		//void seteps2(double ep2){eps2 = ep2;}
		void setPosition(Vector2r p){Position = p;}
		void setOrientation(Rotationr Ori){
		Orientation = Ori;
		//rotation matrices
    rot_mat2local = (Orientation).inverse().toRotationMatrix();//to particle's system
	  rot_mat2global = (Orientation).toRotationMatrix();//to global system
		}

		////////////////////////////////////////////
		Vector2r getSurface(Real phi) const;
		Vector2r getNormal(Real);
		Real Normal2Phi(Vector2r);


    //rotation matrices
    Matrix2r rot_mat2local; //(Orientation).conjugate().toRotationMatrix();//to particle's system
    Matrix2r rot_mat2global; // (Orientation).toRotationMatrix();//to global system


	protected:
		//double rx, ry, rz, eps1, eps2;//the main parameters of a superquadric
		//r_max:the maximum r in the three axes
    double r_max;
    //sign of performed initialization
		bool init;
		Vector2r Position;			//the center position of a superquadric, i.e (x, y, z)
		double Area;
		Real Inertia;			//the pricipal moments of inetia
		Rotationr Orientation;	//orientation
		Vector2r Surface;			//surface equations represented by the sureface parameters Phi1 and Phi2, where Phi1 belongs [-pi,pi], Phi2 in [-pi/2., pi/2.].
		Vector2r CurrentP;			//the local current point (x,y,z) that is used to calculate some parameters.
		Vector2r CurrentPhi;		//the local current phi1,phi2
		Vector2r CurrentNorm;		//the local current n1,n2,n3
		Vector2r ContactNorm;		//the vector of contact direction
public:
        int m_GL_slices;//
        int m_glSolidListNum;//glList number
        int m_glWireListNum;//glList number
		SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(Superellipse,Shape,"Superellipse shape.",
                        ((Vector2r,rxy,,,"length of half-axis in the two driections, x,y, local coordinate system"))
												((Vector2r,ref_rxy,,,"reference length of half-axis in the two driections, x,y, local coordinate system"))
                        ((double,rx,,,"length of half-axis in the two driections, x,y,z, local coordinate system"))
                        ((double,ry,,,"length of half-axis in the two driections, x,y,z, local coordinate system"))
                        ((double,eps,,,"parameter eps in the surfce funtion, i.e., (x/a)^(2/eps)+(y/b)^(2/eps)=1"))
                        ((bool,isSphere,false,,"the shape is non-spherical by default. Using this flag for accelerating spherical systems")),
			/*init*/,
			/*ctor*/
			createIndex();
			init = false;
      m_glSolidListNum = -1;
      m_glWireListNum = -1;
			m_GL_slices = 10,
			.def("Initial",&Superellipse::Initial,"Initialization")
			.def("getArea",&Superellipse::getArea,"return particle's area")
			.def("getrxy",&Superellipse::getrxy,"return particle's rxy")
			.def("geteps",&Superellipse::geteps,"return particle's eps")
			//.def("geteps2",&Superellipse::geteps2,"return particle's eps2")
			.def("getInertia",&Superellipse::getInertia,"return particle's inertia tensor")
			.def("getOri",&Superellipse::getOrientation,"return particle's orientation")
			.def("getCentroid",&Superellipse::getPosition,"return particle's centroid")
		);
		REGISTER_CLASS_INDEX(Superellipse,Shape);
};
REGISTER_SERIALIZABLE(Superellipse);


//***************************************************************************
/*! Collision configuration for Superellipse and something.
 * This is expressed as penetration volume properties: centroid, volume, depth ...
 *
 * Self-contained. */
class SuperellipseGeom: public IGeom{
	public:
		//Real rotAngle;//rotation angle
		Vector2r preNormal;//normal at the previous step
		virtual ~SuperellipseGeom();

		//precompute data for shear evaluation
		void precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector2r& currentNormal, bool isNew, const Vector2r& shift2);
		Vector2r& rotate(Vector2r& shearForce) const;
		//sep_plane is a code storing plane, that previously separated two polyhedras. It is used for faster detection of non-overlap.
		//std::vector<int> sep_plane;
		//bool isShearNew;
	//protected:
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperellipseGeom,IGeom,"Geometry of interaction between 2 :yref:`vector<Superellipse>`, including volumetric characteristics",
			((Real,PenetrationDepth,0.,,"PenetrationDepth"))
			((Real,contactAngle,0.0,,"Normal direction (in local coordinates)"))
			((Vector2r,contactPoint,Vector2r::Zero(),,"Contact point (global coords), center of the PenetrationDepth"))
			((Vector2r,shearInc,Vector2r::Zero(),,"Shear displacement increment in the last step"))
			((Vector2r,relativeVn,Vector2r::Zero(),,"relative velocity in the normal at the contact"))
			((Vector2r,relativeVs,Vector2r::Zero(),,"relative velocity in the tangential at the contact"))
			((Vector2r,normal,Vector2r::Zero(),,"Contact normal"))
			((Vector2r,shift2,Vector2r::Zero(),,""))
			((Vector2r,point1,Vector2r::Zero(),,"point 1 of the closest two points."))
			((Vector2r,point2,Vector2r::Zero(),,"point 2 of the closest two points."))
			((Vector2r,point11,Vector2r::Zero(),,"direction from point 1.FOR TESTING"))
			((Vector2r,point22,Vector2r::Zero(),,"direction from point 2.FOR TESTING"))
			((bool,isShearNew,true,,"")),
			//((Vector2r,orthonormal_axis,Vector2r::Zero(),,"")),
			createIndex();preNormal=Vector2r::Zero();
			//sep_plane.assign(3,0);
		);
		//FUNCTOR2D(Tetra,Tetra);
		REGISTER_CLASS_INDEX(SuperellipseGeom,IGeom);
};
REGISTER_SERIALIZABLE(SuperellipseGeom);

//***************************************************************************
/*! Creates Aabb from Superellipse.
 *
 * Self-contained. */
class Bo1_Superellipse_Aabb: public BoundFunctor{
	public:
		void go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se2r& se2, const Body*);
		FUNCTOR1D(Superellipse);
		SUDODEM_CLASS_BASE_DOC(Bo1_Superellipse_Aabb,BoundFunctor,"Create/update :yref:`Aabb` of a :yref:`Superellipse`");
};
REGISTER_SERIALIZABLE(Bo1_Superellipse_Aabb);

//***************************************************************************
/*! Elastic material */
class SuperellipseMat: public Material{
	public:
		 SuperellipseMat(double N, double S, double F){Kn=N; Ks=S; frictionAngle=F;};
		 double GetStrength(){return strength;};
	virtual ~SuperellipseMat(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperellipseMat,Material,"Elastic material with Coulomb friction.",
		((Real,Kn,1e8,,"Normal stiffness (N/m)."))
		((Real,Ks,1e5,,"Shear stiffness (N/m)."))
		((Real,frictionAngle,.5,,"Contact friction angle (in radians)."))
		((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
		((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
		((bool,IsSplitable,0,,"To be splitted ... or not"))
		((double,strength,100,,"Stress at whis polyhedra of volume 4/3*pi [mm] breaks.")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(SuperellipseMat,Material);
};
REGISTER_SERIALIZABLE(SuperellipseMat);
//material for Hertz-Mindlin model
class SuperellipseMat2: public SuperellipseMat{
	public:
		 SuperellipseMat2(double G_in, double nu_in, double F){G=G_in; nu=nu_in; frictionAngle=F;};
		 //double GetStrength(){return strength;};
	virtual ~SuperellipseMat2(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperellipseMat2,SuperellipseMat,"Elastic material with Coulomb friction.",
		((Real,G,29e9,,"shear modulus (Pa)."))
		((Real,nu,0.15,,"Poisson's ratio.")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(SuperellipseMat2,SuperellipseMat);
};
REGISTER_SERIALIZABLE(SuperellipseMat2);

//use FrictPhys as the basic class
class SuperellipsePhys: public FrictPhys{
	public:
	virtual ~SuperellipsePhys(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperellipsePhys,FrictPhys,"Simple elastic material with friction for volumetric constitutive laws",
		//((Real,kn,0,,"Normal stiffness"))
		//((Vector2r,normalForce,Vector2r::Zero(),,"Normal force after previous step (in global coordinates)."))
		//((Real,ks,0,,"Shear stiffness"))
		//((Vector2r,shearForce,Vector2r::Zero(),,"Shear force after previous step (in global coordinates)."))
		// Contact damping ratio as for linear elastic contact law
		((Vector2r,normalViscous,Vector2r::Zero(),,"Normal viscous component"))
		((Vector2r,shearViscous,Vector2r::Zero(),,"Shear viscous component"))
		((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
		((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$.")),
		//((Real,tangensOfFrictionAngle,0.,,"tangens of angle of internal friction")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(SuperellipsePhys,FrictPhys);
};
REGISTER_SERIALIZABLE(SuperellipsePhys);

//***************************************************************************
#ifdef SUDODEM_OPENGL
	#include<sudodem/pkg/common/GLDrawFunctors.hpp>
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	#include<sudodem/lib/opengl/GLUtils.hpp>
	#include<GL/glu.h>
	#include<sudodem/pkg/dem/Shop.hpp>
	class Gl1_Superellipse : public GlShapeFunctor{
		private:
			// for stripes
			//static int glStripedDiskList;
			//static int glGlutDiskList;
			//Generate GlList for GLUT disk
			void initWireGlList(Superellipse* t);
			//Generate GlList for sliced disks
			void initSolidGlList(Superellipse* t);
			//for regenerating glutDisk list if needed
			static int preSlices;
		public:
			virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
		SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_Superellipse,GlShapeFunctor,"Renders :yref:`Superellipse` object",
			((bool,wire,false,,"Only show wireframe (controlled by ``glutSlices`` and ``glutStacks``."))
			((int,Slices,12,(Attr::noSave | Attr::readonly),"Base number of Superellipse slices"))
		);
		RENDERS(Superellipse);
	};

	REGISTER_SERIALIZABLE(Gl1_Superellipse);

	struct Gl1_SuperellipseGeom: public GlIGeomFunctor{
		RENDERS(SuperellipseGeom);
		void go(const shared_ptr<IGeom>&, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool);
		void draw(const shared_ptr<IGeom>&);
		SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_SuperellipseGeom,GlIGeomFunctor,"Render :yref:`SuperellipseGeom` geometry.",
		);
	};
	REGISTER_SERIALIZABLE(Gl1_SuperellipseGeom);
#endif


//***************************************************************************
class Ip2_SuperellipseMat_SuperellipseMat_SuperellipsePhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
			const shared_ptr<Material>& b2,
			const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(SuperellipseMat,SuperellipseMat);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_SuperellipseMat_SuperellipseMat_SuperellipsePhys,IPhysFunctor,"",
	//((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))//will crash:Abort (core dumped)
	//((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
	);
};
REGISTER_SERIALIZABLE(Ip2_SuperellipseMat_SuperellipseMat_SuperellipsePhys);

//***************************************************************************
/*! Calculate physical response based on penetration configuration given by TTetraGeom. */

class SuperellipseLaw: public LawFunctor{
	OpenMPAccumulator<Real> plasticDissipation;
	virtual bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
	Real elasticEnergy ();
	Real getPlasticDissipation();
	void initPlasticDissipation(Real initVal=0);
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(SuperellipseLaw,LawFunctor,"Calculate physical response of 2 :yref:`vector<Superellipse>` in interaction, based on penetration configuration given by :yref:`SuperellipseGeom`.",
	((Vector2r,shearForce,Vector2r::Zero(),,"Shear force from last step"))
	((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts. This will trace only plastic energy in this law, see O.trackEnergy for a more complete energies tracing"))
	((int,plastDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for plastic dissipation (with O.trackEnergy)"))
	((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
	,,
	.def("elasticEnergy",&SuperellipseLaw::elasticEnergy,"Compute and return the total elastic energy in all \"FrictPhys\" contacts")
	.def("plasticDissipation",&SuperellipseLaw::getPlasticDissipation,"Total energy dissipated in plastic slips at all FrictPhys contacts. Computed only if :yref:`Law2_ScGeom_FrictPhys_CundallStrack::traceEnergy` is true.")
	.def("initPlasticDissipation",&SuperellipseLaw::initPlasticDissipation,"Initialize cummulated plastic dissipation to a value (0 by default).")
	);
	FUNCTOR2D(SuperellipseGeom,SuperellipsePhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SuperellipseLaw);
