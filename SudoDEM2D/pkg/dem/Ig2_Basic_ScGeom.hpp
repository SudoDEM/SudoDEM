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
#include<sudodem/pkg/common/Disk.hpp>
//#include<sudodem/pkg/common/Box.hpp>
//#include<sudodem/pkg/common/Facet.hpp>
#include<sudodem/pkg/common/Wall.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>
//geometry computation between basic shapes:disk-disk, disk-facet,disk-wall,disk-box

class Ig2_Disk_Disk_ScGeom: public IGeomFunctor{
	public:
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c);

		SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Disk_Disk_ScGeom,IGeomFunctor,
		"Create/update a :yref:`ScGeom` instance representing the geometry of a contact point between two :yref:`Disks<Disk>` s.",
		((Real,interactionDetectionFactor,1,,"Enlarge both radii by this factor (if >1), to permit creation of distant interactions.\n\nInteractionGeometry will be computed when interactionDetectionFactor*(rad1+rad2) > distance.\n\n.. note::\n\t This parameter is functionally coupled with :yref:`Bo1_Disk_Aabb::aabbEnlargeFactor`, which will create larger bounding boxes and should be of the same value."))
		((bool,avoidGranularRatcheting,true,,"Define relative velocity so that ratcheting is avoided. It applies for disk-disk contacts. It eventualy also apply for disk-emulating interactions (i.e. convertible into the ScGeom type), if the virtual disk's motion is defined correctly (see e.g. :yref:`Ig2_Disk_ChainedCylinder_CylScGeom`.\n\n"
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
	FUNCTOR2D(Disk,Disk);
	// needed for the dispatcher, even if it is symmetric
	DEFINE_FUNCTOR_ORDER_2D(Disk,Disk);
};
REGISTER_SERIALIZABLE(Ig2_Disk_Disk_ScGeom);

/*
class Ig2_Facet_Disk_ScGeom : public IGeomFunctor
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector2r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector2r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Facet_Disk_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`Facet` and :yref:`Disk`.",
		((Real,shrinkFactor,((void)"no shrinking",0),,"The radius of the inscribed circle of the facet is decreased by the value of the disk's radius multipled by *shrinkFactor*. From the definition of contact point on the surface made of facets, the given surface is not continuous and becomes in effect surface covered with triangular tiles, with gap between the separate tiles equal to the disk's radius multiplied by 2×*shrinkFactor*. If zero, no shrinking is done."))
	);
	DECLARE_LOGGER;
	FUNCTOR2D(Facet,Disk);
	DEFINE_FUNCTOR_ORDER_2D(Facet,Disk);
};

REGISTER_SERIALIZABLE(Ig2_Facet_Disk_ScGeom);

*/

class Ig2_Wall_Disk_ScGeom: public IGeomFunctor{
	public:
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Wall_Disk_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`Wall` and :yref:`Disk`.",
		((bool,noRatch,true,,"Avoid granular ratcheting"))
	);
	FUNCTOR2D(Wall,Disk);
	DEFINE_FUNCTOR_ORDER_2D(Wall,Disk);
};
REGISTER_SERIALIZABLE(Ig2_Wall_Disk_ScGeom);

class Ig2_Fwall_Disk_ScGeom : public IGeomFunctor
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2,const Vector2r& shift2,const bool& force,
					const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2,const Vector2r& shift2,const bool& force,
					const shared_ptr<Interaction>& c);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_Fwall_Disk_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`Fwall` and :yref:`Disk`.",
	);
	DECLARE_LOGGER;
	FUNCTOR2D(Fwall,Disk);
	DEFINE_FUNCTOR_ORDER_2D(Fwall,Disk);
};

REGISTER_SERIALIZABLE(Ig2_Fwall_Disk_ScGeom);


/*
class Ig2_Box_Disk_ScGeom : public IGeomFunctor
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c);

		virtual bool goReverse(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c);

	SUDODEM_CLASS_BASE_DOC(Ig2_Box_Disk_ScGeom,IGeomFunctor,"Create an interaction geometry :yref:`ScGeom` from :yref:`Box` and :yref:`Disk`, representing the box with a projected virtual disk of same radius.")
	FUNCTOR2D(Box,Disk);
	DEFINE_FUNCTOR_ORDER_2D(Box,Disk);
};
REGISTER_SERIALIZABLE(Ig2_Box_Disk_ScGeom);
*/
