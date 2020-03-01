#ifndef GJKPARTICLE_H
#define GJKPARTICLE_H
/*
 * =====================================================================================
 *
 *       Filename:  GJKParticle.h
 *
 *    Description:  GJKParticle have been defined.
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

#include <GL/glut.h>
#include "GJKParticle_shapes.h"
//**********************************************************************************
class GJKParticle: public Shape{
	public:
		//GJKParticle(std::vector<Vector3r> vertices,int shapeType, SD_Scalar radius,SD_Scalar height);
		GJKParticle(std::vector<Vector3r> vertices,int shapeType, SD_Scalar radius,SD_Scalar height,SD_Scalar margin):m_vertices(vertices),m_shapeType(shapeType),m_radius(radius),m_height(height),m_margin(margin)
		{
			createIndex();
			m_init = false;//this is very important
			//initialize to calculate some basic geometric parameters.
			Initial();

		}
		GJKParticle(std::vector<Vector3r> vertices,int shapeType, SD_Scalar radius,SD_Scalar height):m_vertices(vertices),m_shapeType(shapeType),m_radius(radius),m_height(height),m_margin(0.0)
		{
			createIndex();
			//m_shape = NULL;
			m_init = false;//this is very important
			//initialize to calculate some basic geometric parameters.
			Initial();
		}
		virtual~GJKParticle();
		//~GJKParticle(){if (m_shape!=NULL) delete m_shape;}
		//
		//for object

		void setMargin(SD_Scalar margin){m_margin = margin;}
		void setScaling(const Vector3r& scaling){m_xform.scale(scaling);}
		void setPosition(const Vector3r& pos){ 	m_position = pos;m_xform.setOrigin(pos);}
		//void setOrientation(Quaternionr Ori){	m_orientation = Ori;m_xform.setRotation(Ori);}
		void setOrientation(Quaternionr Ori){
		Ori.normalize();m_orientation = Ori;
		//rotation matrices
		rot_mat2local = (m_orientation).conjugate().toRotationMatrix();//to particle's system
		rot_mat2global = (m_orientation).toRotationMatrix();//to global system
		}
		void setSe3(const Se3r& se3){m_xform.setOrigin(se3.position);m_xform.setRotation(se3.orientation);}
		void setMatrix(const float *m){m_xform.setValue(m);assert(m_xform.getBasis().determinant() != SD_Scalar(0.0));}
		void setMatrix(const double *m){m_xform.setValue(m);assert(m_xform.getBasis().determinant() != SD_Scalar(0.0));}
		void getMatrix(float *m) const{m_xform.getValue(m);}
		void getMatrix(double *m) const{m_xform.getValue(m);}

		SD_Scalar supportH(const Vector3r& v) const;
		//virtual Vector3r support(const Vector3r& v) const {	return SUDODEM_PTR_CAST<DT_Convex>(m_shape)->support(v); }
		Vector3r support(const Vector3r& v) const;
		//friend bool intersect(const GJKParticle*, const GJKParticle*, Vector3r& v);

		//friend bool common_point(const GJKParticle*, const GJKParticle*, Vector3r&,  Vector3r&, Vector3r&);

		//friend bool penetration_depth(const GJKParticle*, const GJKParticle*,  Vector3r&, Vector3r&, Vector3r&);

		//friend SD_Scalar closest_points(const GJKParticle*, const GJKParticle*, Vector3r&, Vector3r&);
		//friend bool penetration_depth(const GJKParticle *,const Matrix3r&,  Vector3r&,const GJKParticle*,const Matrix3r&, Vector3r&, Vector3r&, Vector3r&, Vector3r&);

		//for GJK end


		void Initial();				//calculate some basic geometric parameters.
		bool IsInitialized(){return m_init;}
		int getGLslices(){return m_GL_slices;}//prepare data for visualization
		/////////get
	  SD_Scalar getMargin(){return m_margin;}
		Vector3r getPosition(){return m_position;}
		Vector3r getCentroid(){return m_centroid;}//centroid at the global coordinate system when the particle is initialized.
		int getShapeType() const {return m_shapeType;}
		T_MultiIndexBuf getFacetIndexBuf(){return	m_facetIndexBuf;}
		double getVolume(){Initial();return m_volume;}
		double getRadius(){Initial();return m_radius;}
		Vector3r getInertia(){Initial();return m_inertia;}
		Quaternionr getOrientation(){return m_orientation;}

		//void setPosition(Vector3r p){m_position = p;}
		//void setOrientation(Quaternionr Ori){
		//Ori.normalize();m_orientation = Ori;}
		///////////////here inserted
		void ConvexHull(std::vector<Vector3r>);
		//void randomGeometry(int num_vertices);
		void Polyhedron_volume_centroid();
		Matrix3r TetraInertiaTensor(Vector3r av,Vector3r bv,Vector3r cv,Vector3r dv);
		void Polyhedron_inertia();

		////////////////////////
		//rotation matrices
		Matrix3r rot_mat2local; //(Orientation).conjugate().toRotationMatrix();//to particle's system
		Matrix3r rot_mat2global; // (Orientation).toRotationMatrix();//to global system

private:
		T_MultiIndexBuf m_facetIndexBuf;
		Vector3r position;
  	shared_ptr<DT_Convex> m_shape;//FIXME:m_shape was reset unexpectedly and a crash was throwed out. This occurs after AABB tests.
		//another shared_ptr has been adopted (when calling initial()) to clone m_shape to keep the pointer undestroyed. It seems to work fine!
		shared_ptr<DT_Convex> m_shape2;
		//for object
		//SD_Scalar          m_margin;
		MT_Transform m_xform;

		mutable GLuint         m_displayList;
protected:
		bool m_init;

		Vector3r m_position;			//the center position of a superquadric, i.e (x, y, z)
		double m_volume;
		Vector3r m_centroid;		//centriod of particle
		Vector3r m_inertia;			//the pricipal moments of inetia
		Quaternionr m_orientation;	//orientation

public:
		mutable int testflag;
    int m_GL_slices;//
    int m_glSolidListNum;//glList number
    int m_glWireListNum;//glList number
    //triangulation of facets for GLview, indexes of facet points in m_vertices2
		vector<int> facetTri;
        std::vector<Vector3r> m_facetVertices;
		SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(GJKParticle,Shape,"GJKParticle shape.",
            ((int,m_shapeType,0,,"particle shape:sphere(0), polyhedron(1), cone(2), cylinder(3), and cube(4)"))
            ((std::vector<Vector3r>, m_vertices,,,"Input vertices for convex hull in global coordinate system."))
            //((shared_ptr<DT_Shape>,m_shape,,,"pointer of the shape"))
						((std::vector<Vector3r>, m_vertices2,,,"vertices in local coordinate system, and also used to GLview."))
            //((Vector3r,rxyz,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))
            ((SD_Scalar,m_radius,0.0,,"radius for sphere, cone and cylinder."))
						((SD_Scalar,m_height,0.0,,"radius for sphere, cone and cylinder."))
						((SD_Scalar,m_margin,0.0,,"margin used for hybrid-penetration computation."))
			      ((bool,isSphere,false,,"the shape is non-spherical by default. Using this flag for accelerating spherical systems"))
            /*((Vector2r,eps,,,"length of half-axis in the three driections, x,y,z, local coordinate system"))*/,
			/*init*/,
			/*ctor*/
			createIndex();
			m_init = false;
      m_glSolidListNum = -1;
      m_glWireListNum = -1;testflag=0;
			m_GL_slices = 10,
			.def("Initial",&GJKParticle::Initial,"Initialization")
			.def("getVolume",&GJKParticle::getVolume,"return particle's volume")
			//.def("geteps2",&GJKParticle::geteps2,"return particle's eps2")
			.def("getInertia",&GJKParticle::getInertia,"return particle's inertia tensor")
			.def("getOri",&GJKParticle::getOrientation,"return particle's orientation")
			.def("getCentroid",&GJKParticle::getPosition,"return particle's centroid")
		);
		REGISTER_CLASS_INDEX(GJKParticle,Shape);
};
REGISTER_SERIALIZABLE(GJKParticle);






//***************************************************************************
/*! Collision configuration for GJKParticle and something.
 * This is expressed as penetration volume properties: centroid, volume, depth ...
 *
 * Self-contained. */
class GJKParticleGeom: public IGeom{
	public:
		virtual ~GJKParticleGeom();

		//precompute data for shear evaluation
		void precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector3r& currentNormal, bool isNew, const Vector3r& shift2);
		Vector3r& rotate(Vector3r& shearForce) const;
		//sep_plane is a code storing plane, that previously separated two polyhedras. It is used for faster detection of non-overlap.
		//std::vector<int> sep_plane;
		//bool isShearNew;
	protected:
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(GJKParticleGeom,IGeom,"Geometry of interaction between 2 :yref:`vector<GJKParticle>`, including volumetric characteristics",

			((Real,PenetrationDepth,0.,,"PenetrationDepth"))
			((Vector2r,contactAngle,Vector2r::Zero(),,"Normal direction (in local coordinates)"))
			((Vector3r,contactSA,Vector3r::Zero(),,"separating axis vector"))
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
		REGISTER_CLASS_INDEX(GJKParticleGeom,IGeom);
};
REGISTER_SERIALIZABLE(GJKParticleGeom);

//***************************************************************************
/*! Creates Aabb from GJKParticle.
 *
 * Self-contained. */
class Bo1_GJKParticle_Aabb: public BoundFunctor{
	public:
		void go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se3r& se3, const Body*);
		FUNCTOR1D(GJKParticle);
		SUDODEM_CLASS_BASE_DOC(Bo1_GJKParticle_Aabb,BoundFunctor,"Create/update :yref:`Aabb` of a :yref:`GJKParticle`");
};
REGISTER_SERIALIZABLE(Bo1_GJKParticle_Aabb);

//***************************************************************************
/*! Elastic material */
class GJKParticleMat: public Material{
	public:
		 GJKParticleMat(double N, double S, double F){Kn=N; Ks=S; frictionAngle=F;};
		 double GetStrength(){return strength;};
	virtual ~GJKParticleMat(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(GJKParticleMat,Material,"Elastic material with Coulomb friction.",
		((Real,Kn,1e8,,"Normal stiffness (N/m)."))
		((Real,Ks,1e5,,"Shear stiffness (N/m)."))
		((Real,frictionAngle,.5,,"Contact friction angle (in radians)."))
		((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
		((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
		((bool,IsSplitable,0,,"To be splitted ... or not"))
		((double,strength,100,,"Stress at whis polyhedra of volume 4/3*pi [mm] breaks.")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(GJKParticleMat,Material);
};
REGISTER_SERIALIZABLE(GJKParticleMat);

//***************************************************************************
//class GJKParticlePhys: public IPhys{
//	public:
//	virtual ~GJKParticlePhys(){};
//	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(GJKParticlePhys,IPhys,"Simple elastic material with friction for volumetric constitutive laws",
//		((Real,kn,0,,"Normal stiffness"))
//		((Vector3r,normalForce,Vector3r::Zero(),,"Normal force after previous step (in global coordinates)."))
//		((Real,ks,0,,"Shear stiffness"))
//		((Vector3r,shearForce,Vector3r::Zero(),,"Shear force after previous step (in global coordinates)."))
//		((Real,tangensOfFrictionAngle,0.,,"tangens of angle of internal friction")),
//		/*ctor*/ createIndex();
//	);
//	REGISTER_CLASS_INDEX(GJKParticlePhys,IPhys);
//};
//REGISTER_SERIALIZABLE(GJKParticlePhys);
//use FrictPhys as the basic class
class GJKParticlePhys: public FrictPhys{
	public:
	virtual ~GJKParticlePhys(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(GJKParticlePhys,FrictPhys,"Simple elastic material with friction for volumetric constitutive laws",
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
	REGISTER_CLASS_INDEX(GJKParticlePhys,FrictPhys);
};
REGISTER_SERIALIZABLE(GJKParticlePhys);

//***************************************************************************
#ifdef SUDODEM_OPENGL
	#include<sudodem/pkg/common/GLDrawFunctors.hpp>
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	#include<sudodem/lib/opengl/GLUtils.hpp>
	#include<GL/glu.h>
	#include<sudodem/pkg/dem/Shop.hpp>

	/*! Draw GJKParticle using OpenGL */
	class Gl1_GJKParticle: public GlShapeFunctor{
		private:
		    //static int glSolidList;
		    //static int glWireList;
            void initSolidGlList(GJKParticle*);
            void initWireGlList(GJKParticle*);
			//std::vector<Vector3r> gl_vertices;
			//T_MultiIndexBuf gl_facetIndexBuf;
			//static int pre_slices;

		 	//static bool initialized;
		public:
			virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
			SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_GJKParticle,GlShapeFunctor,"Renders :yref:`GJKParticle` object",
			((bool,wire,true,,"Only show wireframe"))
			((int,slices,10,,"Base number of slices."))
			);
			RENDERS(GJKParticle);
	};
	REGISTER_SERIALIZABLE(Gl1_GJKParticle);

	struct Gl1_GJKParticleGeom: public GlIGeomFunctor{
		RENDERS(GJKParticleGeom);
		void go(const shared_ptr<IGeom>&, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool);
		void draw(const shared_ptr<IGeom>&);
		SUDODEM_CLASS_BASE_DOC_STATICATTRS(Gl1_GJKParticleGeom,GlIGeomFunctor,"Render :yref:`GJKParticleGeom` geometry.",
		);
	};
	REGISTER_SERIALIZABLE(Gl1_GJKParticleGeom);
#endif


//***************************************************************************
class Ip2_GJKParticleMat_GJKParticleMat_GJKParticlePhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
			const shared_ptr<Material>& b2,
			const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(GJKParticleMat,GJKParticleMat);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_GJKParticleMat_GJKParticleMat_GJKParticlePhys,IPhysFunctor,"",
	//((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))//will crash:Abort (core dumped)
	//((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
	);
};
REGISTER_SERIALIZABLE(Ip2_GJKParticleMat_GJKParticleMat_GJKParticlePhys);

//***************************************************************************
/*! Calculate physical response based on penetration configuration given by TTetraGeom. */

class GJKParticleLaw: public LawFunctor{
	OpenMPAccumulator<Real> plasticDissipation;
	virtual bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
	Real elasticEnergy ();
	Real getPlasticDissipation();
	void initPlasticDissipation(Real initVal=0);
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(GJKParticleLaw,LawFunctor,"Calculate physical response of 2 :yref:`vector<GJKParticle>` in interaction, based on penetration configuration given by :yref:`GJKParticleGeom`.",
	//((Vector3r,shearForce,Vector3r::Zero(),,"Shear force from last step"))
	((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts. This will trace only plastic energy in this law, see O.trackEnergy for a more complete energies tracing"))
	((int,plastDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for plastic dissipation (with O.trackEnergy)"))
	((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
	,,
	.def("elasticEnergy",&GJKParticleLaw::elasticEnergy,"Compute and return the total elastic energy in all \"FrictPhys\" contacts")
	.def("plasticDissipation",&GJKParticleLaw::getPlasticDissipation,"Total energy dissipated in plastic slips at all FrictPhys contacts. Computed only if :yref:`Law2_ScGeom_FrictPhys_CundallStrack::traceEnergy` is true.")
	.def("initPlasticDissipation",&GJKParticleLaw::initPlasticDissipation,"Initialize cummulated plastic dissipation to a value (0 by default).")
	);
	FUNCTOR2D(GJKParticleGeom,GJKParticlePhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(GJKParticleLaw);
//********************************************************************************

#endif
