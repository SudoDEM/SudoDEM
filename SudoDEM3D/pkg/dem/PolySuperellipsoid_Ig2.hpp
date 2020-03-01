#include"PolySuperellipsoid.hpp"
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/Box.hpp>

//***************************************************************************
/*! Create PolySuperellipsoid (collision geometry) from colliding PolySuperellipsoids. */
class Ig2_PolySuperellipsoid_PolySuperellipsoid_PolySuperellipsoidGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_PolySuperellipsoid_PolySuperellipsoid_PolySuperellipsoidGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(PolySuperellipsoid,PolySuperellipsoid);
		DEFINE_FUNCTOR_ORDER_2D(PolySuperellipsoid,PolySuperellipsoid);
		SUDODEM_CLASS_BASE_DOC(Ig2_PolySuperellipsoid_PolySuperellipsoid_PolySuperellipsoidGeom,IGeomFunctor,"Create/update geometry of collision between 2 PolySuperellipsoids");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_PolySuperellipsoid_PolySuperellipsoid_PolySuperellipsoidGeom);
//***************************************************************************
/*! Create PolySuperellipsoid (collision geometry) from colliding Wall & PolySuperellipsoid. */
class Ig2_Wall_PolySuperellipsoid_PolySuperellipsoidGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Wall_PolySuperellipsoid_PolySuperellipsoidGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Wall,PolySuperellipsoid);
		DEFINE_FUNCTOR_ORDER_2D(Wall,PolySuperellipsoid);
		SUDODEM_CLASS_BASE_DOC(Ig2_Wall_PolySuperellipsoid_PolySuperellipsoidGeom,IGeomFunctor,"Create/update geometry of collision between Wall and PolySuperellipsoid");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Wall_PolySuperellipsoid_PolySuperellipsoidGeom);
//***************************************************************************
/*! Create PolySuperellipsoid (collision geometry) from colliding Facet & PolySuperellipsoid. */
class Ig2_Facet_PolySuperellipsoid_PolySuperellipsoidGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Facet_PolySuperellipsoid_PolySuperellipsoidGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Facet,PolySuperellipsoid);
		DEFINE_FUNCTOR_ORDER_2D(Facet,PolySuperellipsoid);
		SUDODEM_CLASS_BASE_DOC(Ig2_Facet_PolySuperellipsoid_PolySuperellipsoidGeom,IGeomFunctor,"Create/update geometry of collision between Facet and PolySuperellipsoid");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Facet_PolySuperellipsoid_PolySuperellipsoidGeom);
