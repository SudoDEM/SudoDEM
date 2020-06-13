/*************************************************************************
*  Copyright (C) 2020 by Sway Zhao                                       *
*  zhswee@gmail.com                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifndef TRIELEMENT_IG2_H
#define TRIELEMENT_IG2_H

#include "TriElement.hpp"
#include <sudodem/pkg/common/Sphere.hpp>
#include <sudodem/pkg/dem/PolySuperellipsoid.hpp>
#include <sudodem/pkg/dem/GJKParticle.hpp>

#include <sudodem/pkg/dem/FrictPhys.hpp>
#include <sudodem/pkg/dem/RollingResistanceLaw.hpp>

class Ig2_TriElement_Sphere_ScGeom : public IGeomFunctor
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
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_TriElement_Sphere_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`TriElement` and :yref:`Sphere`.",
		((Real,shrinkFactor,((void)"no shrinking",0),,"The radius of the inscribed circle of the facet is decreased by the value of the sphere's radius multipled by *shrinkFactor*. From the definition of contact point on the surface made of facets, the given surface is not continuous and becomes in effect surface covered with triangular tiles, with gap between the separate tiles equal to the sphere's radius multiplied by 2×*shrinkFactor*. If zero, no shrinking is done."))
	);
	DECLARE_LOGGER;
	FUNCTOR2D(TriElement,Sphere);
	DEFINE_FUNCTOR_ORDER_2D(TriElement,Sphere);
};

REGISTER_SERIALIZABLE(Ig2_TriElement_Sphere_ScGeom);

class Ig2_TriElement_PolySuperellipsoid_ScGeom : public IGeomFunctor
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);
		/*virtual bool goReverse(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);*/
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_TriElement_PolySuperellipsoid_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`TriElement` and :yref:`PolySuperellipsoid`.",
		((Real,shrinkFactor,((void)"no shrinking",0),,"The radius of the inscribed circle of the facet is decreased by the value of the sphere's radius multipled by *shrinkFactor*. From the definition of contact point on the surface made of facets, the given surface is not continuous and becomes in effect surface covered with triangular tiles, with gap between the separate tiles equal to the sphere's radius multiplied by 2×*shrinkFactor*. If zero, no shrinking is done."))
	);
	DECLARE_LOGGER;
	FUNCTOR2D(TriElement,PolySuperellipsoid);
	DEFINE_FUNCTOR_ORDER_2D(TriElement,PolySuperellipsoid);
};
REGISTER_SERIALIZABLE(Ig2_TriElement_PolySuperellipsoid_ScGeom);

class Ip2_RolFrictMat_PolySuperellipsoidMat_FrictPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
			const shared_ptr<Material>& b2,
			const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(RolFrictMat,PolySuperellipsoidMat);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_RolFrictMat_PolySuperellipsoidMat_FrictPhys,IPhysFunctor,"Create a :yref:`FrictPhys` from two :yref:`FrictMats<FrictMat>`. The compliance of one sphere under point load is defined here as $1/(E.D)$, with $E$ the stiffness of the sphere and $D$ its diameter. The compliance of the contact itself will be the sum of compliances from each sphere, i.e. $1/(E_1.D_1)+1/(E_2.D_2)$ in the general case, or $2/(E.D)$ in the special case of equal sizes and equal stiffness. Note that summing compliances corresponds to an harmonic average of stiffnesss (as in e.g. [Scholtes2009a]_), which is how kn is actually computed in the :yref:`Ip2_FrictMat_FrictMat_FrictPhys` functor:\n\n $k_n = \\frac{E_1D_1*E_2D_2}{E_1D_1+E_2D_2}=\\frac{k_1*k_2}{k_1+k_2}$, with $k_i=E_iD_i$.\n\n The shear stiffness ks of one sphere is defined via the material parameter :yref:`ElastMat::poisson`, as ks=poisson*kn, and the resulting shear stiffness of the interaction will be also an harmonic average. In the case of a contact between a :yref:`ViscElMat` and a :yref:`FrictMat`, be sure to set :yref:`FrictMat::young` and :yref:`FrictMat::poisson`, otherwise the default value will be used.",
		((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's friction angle. If ``None``, minimum value is used."))
	);
};
REGISTER_SERIALIZABLE(Ip2_RolFrictMat_PolySuperellipsoidMat_FrictPhys);

//GJKParticle and TriElement
class Ig2_TriElement_GJKParticle_ScGeom : public IGeomFunctor
{
	public :
		virtual bool go(const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);
		/*virtual bool goReverse(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);*/
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ig2_TriElement_GJKParticle_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`TriElement` and :yref:`PolySuperellipsoid`.",
		((Real,shrinkFactor,((void)"no shrinking",0),,"The radius of the inscribed circle of the facet is decreased by the value of the sphere's radius multipled by *shrinkFactor*. From the definition of contact point on the surface made of facets, the given surface is not continuous and becomes in effect surface covered with triangular tiles, with gap between the separate tiles equal to the sphere's radius multiplied by 2×*shrinkFactor*. If zero, no shrinking is done."))
	);
	DECLARE_LOGGER;
	FUNCTOR2D(TriElement,GJKParticle);
	DEFINE_FUNCTOR_ORDER_2D(TriElement,GJKParticle);
};
REGISTER_SERIALIZABLE(Ig2_TriElement_GJKParticle_ScGeom);

class Ip2_RolFrictMat_GJKParticleMat_FrictPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
			const shared_ptr<Material>& b2,
			const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(RolFrictMat,GJKParticleMat);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_RolFrictMat_GJKParticleMat_FrictPhys,IPhysFunctor,"Create a :yref:`FrictPhys` from two :yref:`FrictMats<FrictMat>`. The compliance of one sphere under point load is defined here as $1/(E.D)$, with $E$ the stiffness of the sphere and $D$ its diameter. The compliance of the contact itself will be the sum of compliances from each sphere, i.e. $1/(E_1.D_1)+1/(E_2.D_2)$ in the general case, or $2/(E.D)$ in the special case of equal sizes and equal stiffness. Note that summing compliances corresponds to an harmonic average of stiffnesss (as in e.g. [Scholtes2009a]_), which is how kn is actually computed in the :yref:`Ip2_FrictMat_FrictMat_FrictPhys` functor:\n\n $k_n = \\frac{E_1D_1*E_2D_2}{E_1D_1+E_2D_2}=\\frac{k_1*k_2}{k_1+k_2}$, with $k_i=E_iD_i$.\n\n The shear stiffness ks of one sphere is defined via the material parameter :yref:`ElastMat::poisson`, as ks=poisson*kn, and the resulting shear stiffness of the interaction will be also an harmonic average. In the case of a contact between a :yref:`ViscElMat` and a :yref:`FrictMat`, be sure to set :yref:`FrictMat::young` and :yref:`FrictMat::poisson`, otherwise the default value will be used.",
		((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's friction angle. If ``None``, minimum value is used."))
	);
};
REGISTER_SERIALIZABLE(Ip2_RolFrictMat_GJKParticleMat_FrictPhys);
#endif //TRIELEMENT_IG2_H