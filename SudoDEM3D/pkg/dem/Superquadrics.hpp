#pragma once
/*
 * =====================================================================================
 *
 *       Filename:  superquadrics.h
 *
 *    Description:  superquadrics have been defined.
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
class Superquadrics: public Shape{
	public:
		//constructor
		Superquadrics(double x, double y, double z, double ep1, double ep2);//{/*createIndex();*/ rxyz = Vector3r(x,y,z); eps = Vector2r(ep1,ep2); Initial();};
		/*Superquadrics(Vector3r xyz, Vector2r _eps):rxyz(xyz),eps(_eps)
		{
		createIndex();

	        //initialize to calculate some basic geometric parameters.
	        Initial();
	        };*/
		virtual ~Superquadrics();
		void Initial();				//calculate some basic geometric parameters.
		bool IsInitialized(){return init;}
		/////////get
		//double getrx(){return rx;}
		//double getry(){return ry;}
		//double getrz(){return rz;}
		Vector3r getrxyz(){return rxyz;}
		//double geteps1(){return eps1;}
		//double geteps2(){return eps2;}
		Vector2r geteps(){return eps;}
		double getr_max(){return r_max;}
		Vector3r getPosition(){return Position;}
		double getVolume(){Initial();return Volume;}
		Vector3r getInertia(){Initial();return Inertia;}
		Quaternionr getOrientation(){return Orientation;}
        bool isInside(Vector3r p);
        //curvature related
	    Vector2r Curvatures(Vector2r phi,Vector3r &CurvatureDir_max,Vector3r &CurvatureDir_min,bool &curvatureMaxflag);//principal curvatures


		//void setrx(double x){rx = x;}
		//void setry(double y){ry = y;}
		//void setrz(double z){rz = z;}
		//void seteps1(double ep1){eps1 = ep1;}
		//void seteps2(double ep2){eps2 = ep2;}
		void setPosition(Vector3r p){Position = p;}
		void setOrientation(Quaternionr Ori){
		Ori.normalize();Orientation = Ori;
		//rotation matrices
                rot_mat2local = (Orientation).conjugate().toRotationMatrix();//to particle's system
	        rot_mat2global = (Orientation).toRotationMatrix();//to global system
		}

		////////////////////////////////////////////
		Vector3r getSurface(Vector2r phi) const;
		Vector3r getNormal(Vector2r);
		Vector2r Normal2Phi(Vector3r);
		void getCurrentP(Vector2r);      //get the coordinates of current point P.
		//deraivates
		Mat32r Derivate_P2Phi(Vector2r phi);  //3*2 matrix
		Vector3r Derivate_P2Phi_12(Vector2r phi);
		Vector3r Derivate_P2Phi_11(Vector2r phi);
		Vector3r Derivate_P2Phi_22(Vector2r phi);


		Mat23r Derivate_Phi2N(Vector2r phi, Vector3r n);
		Mat32r Derivate_C2Alpha(Vector2r alpha, Matrix3r RotationMat);
		//fuctions for contact detection
		Vector3r P_alpha12(Vector2r p, Matrix3r globalContact, int sign);
		Vector3r P_alpha12(Vector2r p, Vector3r &n, Matrix3r globalContact,Matrix3r RotationMat, int sign);   //surface represented with alpha1 and alpha2. Alpha1 and alpha2 determine the contact direction, i.e., contect vector c.
		Vector3r P_alpha12(Vector2r para,Matrix3r globalContact,Matrix3r RotationMat,int sign);
		Mat32r JacFunc(Vector2r para, Matrix3r globalContact, Matrix3r RotationMat, int sign);
		Mat32r JacFunc(Vector2r p, Matrix3r globalContact, int sign);		//Jacbian matrix
		Vector3r N_alpha12(Vector2r para, Matrix3r globalContact,int sign);
                //rotation matrices
                Matrix3r rot_mat2local; //(Orientation).conjugate().toRotationMatrix();//to particle's system
	        Matrix3r rot_mat2global; // (Orientation).toRotationMatrix();//to global system


	protected:
		//double rx, ry, rz, eps1, eps2;//the main parameters of a superquadric
		//r_max:the maximum r in the three axes
	        double r_max;
	        //sign of performed initialization
		bool init;
		Vector3r Position;			//the center position of a superquadric, i.e (x, y, z)
		double Volume;
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
		SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(Superquadrics,Shape,"Superquadrics shape.",
                        ((Vector3r,rxyz,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
                        ((double,rx,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
                        ((double,ry,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
                        ((double,rz,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
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
			.def("Initial",&Superquadrics::Initial,"Initialization")
			.def("getVolume",&Superquadrics::getVolume,"return particle's volume")
			.def("getrxyz",&Superquadrics::getrxyz,"return particle's rxyz")
			//.def("getry",&Superquadrics::getry,"return particle's ry")
			//.def("getrz",&Superquadrics::getrz,"return particle's rz")
			.def("geteps",&Superquadrics::geteps,"return particle's eps1 and eps2")
			//.def("geteps2",&Superquadrics::geteps2,"return particle's eps2")
			.def("getInertia",&Superquadrics::getInertia,"return particle's inertia tensor")
			.def("getOri",&Superquadrics::getOrientation,"return particle's orientation")
			.def("getCentroid",&Superquadrics::getPosition,"return particle's centroid")
		);
		REGISTER_CLASS_INDEX(Superquadrics,Shape);
};
REGISTER_SERIALIZABLE(Superquadrics);


//***************************************************************************
/*! Collision configuration for Superquadrics and something.
 * This is expressed as penetration volume properties: centroid, volume, depth ...
 *
 * Self-contained. */
class SuperquadricsGeom: public IGeom{
	public:
		virtual ~SuperquadricsGeom();

		//precompute data for shear evaluation
		void precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector3r& currentNormal, bool isNew, const Vector3r& shift2);
		Vector3r& rotate(Vector3r& shearForce) const;
		//sep_plane is a code storing plane, that previously separated two polyhedras. It is used for faster detection of non-overlap.
		//std::vector<int> sep_plane;
		//bool isShearNew;
	//protected:
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperquadricsGeom,IGeom,"Geometry of interaction between 2 :yref:`vector<Superquadrics>`, including volumetric characteristics",

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
		REGISTER_CLASS_INDEX(SuperquadricsGeom,IGeom);
};
REGISTER_SERIALIZABLE(SuperquadricsGeom);

//SuperquadricsGeom2 is for Hertz-Mindlin model.
class SuperquadricsGeom2: public SuperquadricsGeom{
	public:
		virtual ~SuperquadricsGeom2();
        void HertzMindlin(Vector2r curvature1, Vector2r curvature2, double theta);
	protected:
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperquadricsGeom2,SuperquadricsGeom,"Geometry of interaction between 2 :yref:`vector<Superquadrics>`, including volumetric characteristics",

			((Vector2r,curvatures1,Vector2r::Zero(),,"Principal curvatures of particle 1"))
            ((Vector2r,curvatures2,Vector2r::Zero(),,"Principal curvatures of particle 2"))
            ((Real,curvatureAngle,0.,,"angle between the two maximum principal curvatures"))
            ((Real,alpha,0.,,"semi-length ratio of the contact ellipse,alpha>=1"))
            ((Real,K,0.,,"the complete ellipitc integral of the first kind with argument of e"))
            ((Real,E,0.,,"the complete ellipitc integral of the second kind with argument of e"))
            ((Real,curvatures_sum,0.,,"sum of relative curvatures, A+B"))
			((bool,isSphere,false,,"are the two touching bodies spheres")),
			createIndex();
			//sep_plane.assign(3,0);
		);
		//FUNCTOR2D(Tetra,Tetra);
		REGISTER_CLASS_INDEX(SuperquadricsGeom2,SuperquadricsGeom);
};
REGISTER_SERIALIZABLE(SuperquadricsGeom2);
//***************************************************************************
/*! Creates Aabb from Superquadrics.
 *
 * Self-contained. */
class Bo1_Superquadrics_Aabb: public BoundFunctor{
	public:
		void go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se3r& se3, const Body*);
		FUNCTOR1D(Superquadrics);
		SUDODEM_CLASS_BASE_DOC(Bo1_Superquadrics_Aabb,BoundFunctor,"Create/update :yref:`Aabb` of a :yref:`Superquadrics`");
};
REGISTER_SERIALIZABLE(Bo1_Superquadrics_Aabb);

//***************************************************************************
/*! Elastic material */
class SuperquadricsMat: public Material{
	public:
		 SuperquadricsMat(double N, double S, double F){Kn=N; Ks=S; frictionAngle=F;};
		 double GetStrength(){return strength;};
	virtual ~SuperquadricsMat(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperquadricsMat,Material,"Elastic material with Coulomb friction.",
		((Real,Kn,1e8,,"Normal stiffness (N/m)."))
		((Real,Ks,1e5,,"Shear stiffness (N/m)."))
		((Real,frictionAngle,.5,,"Contact friction angle (in radians)."))
		((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
		((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
		((bool,IsSplitable,0,,"To be splitted ... or not"))
		((double,strength,100,,"Stress at whis polyhedra of volume 4/3*pi [mm] breaks.")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(SuperquadricsMat,Material);
};
REGISTER_SERIALIZABLE(SuperquadricsMat);
//material for Hertz-Mindlin model
class SuperquadricsMat2: public SuperquadricsMat{
	public:
		 SuperquadricsMat2(double G_in, double nu_in, double F){G=G_in; nu=nu_in; frictionAngle=F;};
		 //double GetStrength(){return strength;};
	virtual ~SuperquadricsMat2(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperquadricsMat2,SuperquadricsMat,"Elastic material with Coulomb friction.",
		((Real,G,29e9,,"shear modulus (Pa)."))
		((Real,nu,0.15,,"Poisson's ratio.")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(SuperquadricsMat2,SuperquadricsMat);
};
REGISTER_SERIALIZABLE(SuperquadricsMat2);


//***************************************************************************

//class SuperquadricsPhys: public IPhys{
//	public:
//	virtual ~SuperquadricsPhys(){};
//	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperquadricsPhys,IPhys,"Simple elastic material with friction for volumetric constitutive laws",
//		((Real,kn,0,,"Normal stiffness"))
//		((Vector3r,normalForce,Vector3r::Zero(),,"Normal force after previous step (in global coordinates)."))
//		((Real,ks,0,,"Shear stiffness"))
//		((Vector3r,shearForce,Vector3r::Zero(),,"Shear force after previous step (in global coordinates)."))
//		// Contact damping ratio as for linear elastic contact law
//		((Vector3r,normalViscous,Vector3r::Zero(),,"Normal viscous component"))
//		((Vector3r,shearViscous,Vector3r::Zero(),,"Shear viscous component"))
//		((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
//		((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
//		((Real,tangensOfFrictionAngle,0.,,"tangens of angle of internal friction")),
//		/*ctor*/ createIndex();
//	);
//	REGISTER_CLASS_INDEX(SuperquadricsPhys,IPhys);
//};
//REGISTER_SERIALIZABLE(SuperquadricsPhys);
//use FrictPhys as the basic class
class SuperquadricsPhys: public FrictPhys{
	public:
	virtual ~SuperquadricsPhys(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperquadricsPhys,FrictPhys,"Simple elastic material with friction for volumetric constitutive laws",
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
	REGISTER_CLASS_INDEX(SuperquadricsPhys,FrictPhys);
};
REGISTER_SERIALIZABLE(SuperquadricsPhys);

//////////////////////Herzian contact law
class SuperquadricsHertzMindlinPhys: public SuperquadricsPhys{
	public:
	virtual ~SuperquadricsHertzMindlinPhys(){};

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SuperquadricsHertzMindlinPhys,SuperquadricsPhys,"Simple elastic material with friction for volumetric constitutive laws",
		((Real,nuab,0,,"average nu of two particles"))
		((Real,Gab,0,,"Shear modulus"))
		((Real,mu1,0.,,"parameter"))
		((Real,mu2,0.,,"parameter")),
		//((bool,isNew,true,,"parameter")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(SuperquadricsHertzMindlinPhys,SuperquadricsPhys);
};
REGISTER_SERIALIZABLE(SuperquadricsHertzMindlinPhys);

//***************************************************************************
#ifdef SUDODEM_OPENGL
	#include<sudodem/pkg/common/GLDrawFunctors.hpp>
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	#include<sudodem/lib/opengl/GLUtils.hpp>
	#include<GL/glu.h>
	#include<sudodem/pkg/dem/Shop.hpp>

	/*! Draw Superquadrics using OpenGL */
	class Gl1_Superquadrics: public GlShapeFunctor{
		private:
			//static vector<Vector3r> vertices;
            static int glSolidList;
		    static int glWireList;
            void initSolidGlList(Superquadrics*);
            void initWireGlList(Superquadrics*);
            static int pre_slices;
			void drawSlice(Vector3r &p1,Vector3r &p2,Vector3r &p3,Vector3r &p4);

		public:
			virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
			SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_Superquadrics,GlShapeFunctor,"Renders :yref:`Superquadrics` object",
			((bool,wire,true,,"Only show wireframe"))
			((int,slices,5,,"Base number of slices."))
			);
			RENDERS(Superquadrics);
	};
	REGISTER_SERIALIZABLE(Gl1_Superquadrics);

	struct Gl1_SuperquadricsGeom: public GlIGeomFunctor{
		RENDERS(SuperquadricsGeom);
		void go(const shared_ptr<IGeom>&, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool);
		void draw(const shared_ptr<IGeom>&);
		SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_SuperquadricsGeom,GlIGeomFunctor,"Render :yref:`SuperquadricsGeom` geometry.",
		);
	};
	REGISTER_SERIALIZABLE(Gl1_SuperquadricsGeom);
#endif


//***************************************************************************
class Ip2_SuperquadricsMat_SuperquadricsMat_SuperquadricsPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
			const shared_ptr<Material>& b2,
			const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(SuperquadricsMat,SuperquadricsMat);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_SuperquadricsMat_SuperquadricsMat_SuperquadricsPhys,IPhysFunctor,"",
	//((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))//will crash:Abort (core dumped)
	//((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
	);
};
REGISTER_SERIALIZABLE(Ip2_SuperquadricsMat_SuperquadricsMat_SuperquadricsPhys);
//***************************************************************************
//Hertz-Mindlin model
class Ip2_SuperquadricsMat2_SuperquadricsMat2_HertzMindlinPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
			const shared_ptr<Material>& b2,
			const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(SuperquadricsMat2,SuperquadricsMat2);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_SuperquadricsMat2_SuperquadricsMat2_HertzMindlinPhys,IPhysFunctor,"",
	//((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))//will crash:Abort (core dumped)
	//((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
	);
};
REGISTER_SERIALIZABLE(Ip2_SuperquadricsMat2_SuperquadricsMat2_HertzMindlinPhys);

//***************************************************************************
/*! Calculate physical response based on penetration configuration given by TTetraGeom. */

class SuperquadricsLaw: public LawFunctor{
	OpenMPAccumulator<Real> plasticDissipation;
	virtual bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
	Real elasticEnergy ();
	Real getPlasticDissipation();
	void initPlasticDissipation(Real initVal=0);
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(SuperquadricsLaw,LawFunctor,"Calculate physical response of 2 :yref:`vector<Superquadrics>` in interaction, based on penetration configuration given by :yref:`SuperquadricsGeom`.",
	((Vector3r,shearForce,Vector3r::Zero(),,"Shear force from last step"))
	((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts. This will trace only plastic energy in this law, see O.trackEnergy for a more complete energies tracing"))
	((int,plastDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for plastic dissipation (with O.trackEnergy)"))
	((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
	,,
	.def("elasticEnergy",&SuperquadricsLaw::elasticEnergy,"Compute and return the total elastic energy in all \"FrictPhys\" contacts")
	.def("plasticDissipation",&SuperquadricsLaw::getPlasticDissipation,"Total energy dissipated in plastic slips at all FrictPhys contacts. Computed only if :yref:`Law2_ScGeom_FrictPhys_CundallStrack::traceEnergy` is true.")
	.def("initPlasticDissipation",&SuperquadricsLaw::initPlasticDissipation,"Initialize cummulated plastic dissipation to a value (0 by default).")
	);
	FUNCTOR2D(SuperquadricsGeom,SuperquadricsPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SuperquadricsLaw);

class SuperquadricsLaw2: public LawFunctor{

	virtual bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(SuperquadricsLaw2,LawFunctor,"Calculate physical response of 2 :yref:`vector<Superquadrics>` in interaction, based on penetration configuration given by :yref:`SuperquadricsGeom`.",
	((Vector3r,shearForce,Vector3r::Zero(),,"Shear force from last step"))
	,,

	);
	FUNCTOR2D(SuperquadricsGeom2,SuperquadricsHertzMindlinPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SuperquadricsLaw2);
/*
class SuperquadricsLaw2: public SuperquadricsLaw{

	virtual bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(SuperquadricsLaw2,LawFunctor,"Calculate physical response of 2 :yref:`vector<Superquadrics>` in interaction, based on penetration configuration given by :yref:`SuperquadricsGeom`.",
	,,

	);
	FUNCTOR2D(SuperquadricsGeom2,SuperquadricsHertzMindlinPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SuperquadricsLaw2);
*/
//********************************************************************************
shared_ptr<Body> NewSuperquadrics(double x, double y, double z, double ep1, double ep2, shared_ptr<Material> mat,bool rotate);
//shared_ptr<Body> NewSuperquadrics(Vector3r rxyz, Vector2r eps, shared_ptr<Material> mat,bool rotate);
