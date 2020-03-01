/*************************************************************************
*  Copyright (C) 2018 by Shiwei Zhao                                     *
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*  Copyright (C) 2006 by Bruno Chareyre                                  *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#include<sudodem/lib/serialization/Serializable.hpp>
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/Box.hpp>
#include<sudodem/pkg/common/Facet.hpp>
#include<sudodem/pkg/common/Wall.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>
//geometry computation between basic shapes:sphere-sphere, sphere-facet,sphere-wall,sphere-box

class Ig2_Sphere_Sphere_ScGeom: public IGeomFunctor{
	public:
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);

		SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Sphere_Sphere_ScGeom,IGeomFunctor,
		"Create/update a :yref:`ScGeom` instance representing the geometry of a contact point between two :yref:`Spheres<Sphere>` s.",
		((Real,interactionDetectionFactor,1,,"Enlarge both radii by this factor (if >1), to permit creation of distant interactions.\n\nInteractionGeometry will be computed when interactionDetectionFactor*(rad1+rad2) > distance.\n\n.. note::\n\t This parameter is functionally coupled with :yref:`Bo1_Sphere_Aabb::aabbEnlargeFactor`, which will create larger bounding boxes and should be of the same value."))
		((bool,avoidGranularRatcheting,true,,"Define relative velocity so that ratcheting is avoided. It applies for sphere-sphere contacts. It eventualy also apply for sphere-emulating interactions (i.e. convertible into the ScGeom type), if the virtual sphere's motion is defined correctly (see e.g. :yref:`Ig2_Sphere_ChainedCylinder_CylScGeom`.\n\n"
		"Short explanation of what we want to avoid :\n\n"
		"Numerical ratcheting is best understood considering a small elastic cycle at a contact between two grains : assuming b1 is fixed, impose this displacement to b2 :\n\n"
  		"#. translation *dx* in the normal direction\n"
		"#. rotation *a*\n"
		"#. translation *-dx* (back to the initial position)\n"
		"#. rotation *-a* (back to the initial orientation)\n\n\n"
		"If the branch vector used to define the relative shear in rotation×branch is not constant (typically if it is defined from the vector center→contactPoint), then the shear displacement at the end of this cycle is not zero: rotations *a* and *-a* are multiplied by branches of different lengths.\n\n"
		"It results in a finite contact force at the end of the cycle even though the positions and orientations are unchanged, in total contradiction with the elastic nature of the problem. It could also be seen as an *inconsistent energy creation or loss*. Given that DEM simulations tend to generate oscillations around equilibrium (damped mass-spring), it can have a significant impact on the evolution of the packings, resulting for instance in slow creep in iterations under constant load.\n\n"
		"The solution adopted here to avoid ratcheting is as proposed by McNamara and co-workers. They analyzed the ratcheting problem in detail - even though they comment on the basis of a cycle that differs from the one shown above. One will find interesting discussions in e.g. [McNamara2008]_, even though solution it suggests is not fully applied here (equations of motion are not incorporating alpha, in contradiction with what is suggested by McNamara et al.).\n\n"
		))
	);
	FUNCTOR2D(Sphere,Sphere);
	// needed for the dispatcher, even if it is symmetric
	DEFINE_FUNCTOR_ORDER_2D(Sphere,Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Sphere_Sphere_ScGeom);

class Ig2_Sphere_Sphere_ScGeom6D: public Ig2_Sphere_Sphere_ScGeom{
	public:
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);

		SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Sphere_Sphere_ScGeom6D,Ig2_Sphere_Sphere_ScGeom,"Create/update a :yref:`ScGeom6D` instance representing the geometry of a contact point between two :yref:`Spheres<Sphere>`, including relative rotations.",
		((bool,updateRotations,true,,"Precompute relative rotations. Turning this false can speed up simulations when rotations are not needed in constitutive laws (e.g. when spheres are compressed without cohesion and moment in early stage of a triaxial test), but is not foolproof. Change this value only if you know what you are doing."))
		((bool,creep,false,,"Substract rotational creep from relative rotation. The rotational creep :yref:`ScGeom6D::twistCreep` is a quaternion and has to be updated inside a constitutive law, see for instance :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`."
		))
	);
	FUNCTOR2D(Sphere,Sphere);
	// needed for the dispatcher, even if it is symmetric
	DEFINE_FUNCTOR_ORDER_2D(Sphere,Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Sphere_Sphere_ScGeom6D);

class Ig2_Facet_Sphere_ScGeom : public IGeomFunctor
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Facet_Sphere_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`Facet` and :yref:`Sphere`.",
		((Real,shrinkFactor,((void)"no shrinking",0),,"The radius of the inscribed circle of the facet is decreased by the value of the sphere's radius multipled by *shrinkFactor*. From the definition of contact point on the surface made of facets, the given surface is not continuous and becomes in effect surface covered with triangular tiles, with gap between the separate tiles equal to the sphere's radius multiplied by 2×*shrinkFactor*. If zero, no shrinking is done."))
	);
	DECLARE_LOGGER;
	FUNCTOR2D(Facet,Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Facet,Sphere);
};

REGISTER_SERIALIZABLE(Ig2_Facet_Sphere_ScGeom);

class Ig2_Facet_Sphere_ScGeom6D : public Ig2_Facet_Sphere_ScGeom
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);
	SUDODEM_CLASS_BASE_DOC(Ig2_Facet_Sphere_ScGeom6D,Ig2_Facet_Sphere_ScGeom,"Create an interaction geometry :yref:`ScGeom6D` from :yref:`Facet` and :yref:`Sphere`, representing the Facet with a projected virtual sphere of same radius.")
	FUNCTOR2D(Facet,Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Facet,Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Facet_Sphere_ScGeom6D);

class Ig2_Wall_Sphere_ScGeom: public IGeomFunctor{
	public:
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Wall_Sphere_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`Wall` and :yref:`Sphere`.",
		((bool,noRatch,true,,"Avoid granular ratcheting"))
	);
	FUNCTOR2D(Wall,Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Wall,Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Wall_Sphere_ScGeom);


class Ig2_Box_Sphere_ScGeom : public IGeomFunctor
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);

		virtual bool goReverse(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);

	SUDODEM_CLASS_BASE_DOC(Ig2_Box_Sphere_ScGeom,IGeomFunctor,"Create an interaction geometry :yref:`ScGeom` from :yref:`Box` and :yref:`Sphere`, representing the box with a projected virtual sphere of same radius.")
	FUNCTOR2D(Box,Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Box,Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Box_Sphere_ScGeom);

class Ig2_Box_Sphere_ScGeom6D : public Ig2_Box_Sphere_ScGeom
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);

		virtual bool goReverse(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);

	SUDODEM_CLASS_BASE_DOC(Ig2_Box_Sphere_ScGeom6D,Ig2_Box_Sphere_ScGeom,"Create an interaction geometry :yref:`ScGeom6D` from :yref:`Box` and :yref:`Sphere`, representing the box with a projected virtual sphere of same radius.")
	FUNCTOR2D(Box,Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Box,Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Box_Sphere_ScGeom6D);
