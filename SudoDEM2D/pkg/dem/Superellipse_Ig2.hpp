#include"Superellipse.hpp"
//#include<sudodem/pkg/common/Sphere.hpp>
//#include<sudodem/pkg/common/Box.hpp>

//***************************************************************************
/*! Create Superellipse (collision geometry) from colliding Superellipses. */
class Ig2_Superellipse_Superellipse_SuperellipseGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Superellipse_Superellipse_SuperellipseGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Superellipse,Superellipse);
		DEFINE_FUNCTOR_ORDER_2D(Superellipse,Superellipse);
		SUDODEM_CLASS_BASE_DOC(Ig2_Superellipse_Superellipse_SuperellipseGeom,IGeomFunctor,"Create/update geometry of collision between 2 Superellipses");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Superellipse_Superellipse_SuperellipseGeom);
//***************************************************************************
/*! Create Superellipse (collision geometry) from colliding Wall & Superellipse. */
class Ig2_Wall_Superellipse_SuperellipseGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Wall_Superellipse_SuperellipseGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Wall,Superellipse);
		DEFINE_FUNCTOR_ORDER_2D(Wall,Superellipse);
		SUDODEM_CLASS_BASE_DOC(Ig2_Wall_Superellipse_SuperellipseGeom,IGeomFunctor,"Create/update geometry of collision between Wall and Superellipse");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Wall_Superellipse_SuperellipseGeom);
//***************************************************************************
class Ig2_Fwall_Superellipse_SuperellipseGeom : public IGeomFunctor
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2,const Vector2r& shift2,const bool& force,
					const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2,const Vector2r& shift2,const bool& force,
					const shared_ptr<Interaction>& c);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Fwall_Superellipse_SuperellipseGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`Fwall` and :yref:`Disk`.",
	);
	DECLARE_LOGGER;
	FUNCTOR2D(Fwall,Superellipse);
	DEFINE_FUNCTOR_ORDER_2D(Fwall,Superellipse);
};

REGISTER_SERIALIZABLE(Ig2_Fwall_Superellipse_SuperellipseGeom);
//***************************************************************************
/*! Create Superellipse (collision geometry) from colliding Facet & Superellipse. */
/*
class Ig2_Facet_Superellipse_SuperellipseGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Facet_Superellipse_SuperellipseGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Facet,Superellipse);
		DEFINE_FUNCTOR_ORDER_2D(Facet,Superellipse);
		SUDODEM_CLASS_BASE_DOC(Ig2_Facet_Superellipse_SuperellipseGeom,IGeomFunctor,"Create/update geometry of collision between Facet and Superellipse");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Facet_Superellipse_SuperellipseGeom);

*/
