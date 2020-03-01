#include"Superquadrics.hpp"
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/Box.hpp>

//***************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Superquadricss. */
class Ig2_Superquadrics_Superquadrics_SuperquadricsGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Superquadrics_Superquadrics_SuperquadricsGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Superquadrics,Superquadrics);
		DEFINE_FUNCTOR_ORDER_2D(Superquadrics,Superquadrics);
		SUDODEM_CLASS_BASE_DOC(Ig2_Superquadrics_Superquadrics_SuperquadricsGeom,IGeomFunctor,"Create/update geometry of collision between 2 Superquadricss");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Superquadrics_Superquadrics_SuperquadricsGeom);
//***************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Superquadricss. */
//This is for Hertz-Mindlin model
class Ig2_Superquadrics_Superquadrics_SuperquadricsGeom2: public IGeomFunctor
{
	public:
		virtual ~Ig2_Superquadrics_Superquadrics_SuperquadricsGeom2(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Superquadrics,Superquadrics);
		DEFINE_FUNCTOR_ORDER_2D(Superquadrics,Superquadrics);
		SUDODEM_CLASS_BASE_DOC(Ig2_Superquadrics_Superquadrics_SuperquadricsGeom2,IGeomFunctor,"Create/update geometry of collision between 2 Superquadricss");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Superquadrics_Superquadrics_SuperquadricsGeom2);

//***************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Wall & Superquadrics. */
class Ig2_Wall_Superquadrics_SuperquadricsGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Wall_Superquadrics_SuperquadricsGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Wall,Superquadrics);
		DEFINE_FUNCTOR_ORDER_2D(Wall,Superquadrics);
		SUDODEM_CLASS_BASE_DOC(Ig2_Wall_Superquadrics_SuperquadricsGeom,IGeomFunctor,"Create/update geometry of collision between Wall and Superquadrics");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Wall_Superquadrics_SuperquadricsGeom);
//***************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Wall & Superquadrics. */
class Ig2_Wall_Superquadrics_SuperquadricsGeom2: public IGeomFunctor
{
	public:
		virtual ~Ig2_Wall_Superquadrics_SuperquadricsGeom2(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Wall,Superquadrics);
		DEFINE_FUNCTOR_ORDER_2D(Wall,Superquadrics);
		SUDODEM_CLASS_BASE_DOC(Ig2_Wall_Superquadrics_SuperquadricsGeom2,IGeomFunctor,"Create/update geometry of collision between Wall and Superquadrics");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Wall_Superquadrics_SuperquadricsGeom2);
//***************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Facet & Superquadrics. */
class Ig2_Facet_Superquadrics_SuperquadricsGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Facet_Superquadrics_SuperquadricsGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Facet,Superquadrics);
		DEFINE_FUNCTOR_ORDER_2D(Facet,Superquadrics);
		SUDODEM_CLASS_BASE_DOC(Ig2_Facet_Superquadrics_SuperquadricsGeom,IGeomFunctor,"Create/update geometry of collision between Facet and Superquadrics");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Facet_Superquadrics_SuperquadricsGeom);

