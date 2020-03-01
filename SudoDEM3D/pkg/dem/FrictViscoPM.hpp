/*************************************************************************
*  Copyright (C) 2014 by Klaus Thoeni                                    *
*  klaus.thoeni@newcastle.edu.au                                         *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

/**
=== OVERVIEW OF FrictViscoPM ===

A particle model for friction and viscous damping in normal direction. The damping coefficient
has to be defined with the material whereas the contact stiffness kn and ks/kn can be defined with the Ip2 functor.

Remarks:
- maybe there is a better way of implementing this without copying from ElasticContactLaw
- maybe we can combine some ideas of this contact law with other contact laws
*/

#pragma once

#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>
#include<sudodem/pkg/common/MatchMaker.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/lib/base/openmp-accu.hpp>

/** This class holds information associated with each body */
class FrictViscoMat: public FrictMat {
	public:
		virtual ~FrictViscoMat();
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(FrictViscoMat,FrictMat,"Material for use with the FrictViscoPM classes",
			((Real,betan,0.,,"Fraction of the viscous damping coefficient in normal direction equal to $\\frac{c_{n}}{C_{n,crit}}$."))
			,
			createIndex();
		);
		DECLARE_LOGGER;
		REGISTER_CLASS_INDEX(FrictViscoMat,FrictMat);
};
REGISTER_SERIALIZABLE(FrictViscoMat);


/** This class holds information associated with each interaction */
class FrictViscoPhys: public FrictPhys {
	public:
		virtual ~FrictViscoPhys();
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(FrictViscoPhys,FrictPhys,"Representation of a single interaction of the FrictViscoPM type, storage for relevant parameters",
			((Real,cn_crit,NaN,,"Normal viscous constant for ctitical damping defined as $\\c_{n}=C_{n,crit}\\beta_n$."))
			((Real,cn,NaN,,"Normal viscous constant defined as $\\c_{n}=c_{n,crit}\\beta_n$."))
			((Vector3r,normalViscous,Vector3r::Zero(),,"Normal viscous component"))
			,
			createIndex();
			,
		);
		DECLARE_LOGGER;
		REGISTER_CLASS_INDEX(FrictViscoPhys,FrictPhys);
};
REGISTER_SERIALIZABLE(FrictViscoPhys);

/** 2d functor creating IPhys (Ip2) taking FrictViscoMat and FrictViscoMat of 2 bodies, returning type FrictViscoPhys */
class Ip2_FrictViscoMat_FrictViscoMat_FrictViscoPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& pp1, const shared_ptr<Material>& pp2, const shared_ptr<Interaction>& interaction);

		FUNCTOR2D(FrictViscoMat,FrictViscoMat);

		SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_FrictViscoMat_FrictViscoMat_FrictViscoPhys,IPhysFunctor,"Converts 2 :yref:`FrictViscoMat` instances to :yref:`FrictViscoPhys` with corresponding parameters. Basically this functor corresponds to :yref:`Ip2_FrictMat_FrictMat_FrictPhys` with the only difference that damping in normal direction can be considered.",
		((shared_ptr<MatchMaker>,kn,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's normal contact stiffnesses. If this value is not given the elastic properties (i.e. young) of the two colliding materials are used to calculate the stiffness."))
		((shared_ptr<MatchMaker>,kRatio,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's shear contact stiffnesses. If this value is not given the elastic properties (i.e. poisson) of the two colliding materials are used to calculate the stiffness."))
		((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's friction angle. If ``None``, minimum value is used."))
		);
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictViscoMat_FrictViscoMat_FrictViscoPhys);

/** 2d functor creating IPhys (Ip2) taking FrictMat and FrictViscoMat of 2 bodies, returning type FrictViscoPhys */
class Ip2_FrictMat_FrictViscoMat_FrictViscoPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& pp1, const shared_ptr<Material>& pp2, const shared_ptr<Interaction>& interaction);

		FUNCTOR2D(FrictMat,FrictViscoMat);
		DEFINE_FUNCTOR_ORDER_2D(FrictMat,FrictViscoMat);

		SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_FrictMat_FrictViscoMat_FrictViscoPhys,IPhysFunctor,"Converts a :yref:`FrictMat` and :yref:`FrictViscoMat` instance to :yref:`FrictViscoPhys` with corresponding parameters. Basically this functor corresponds to :yref:`Ip2_FrictMat_FrictMat_FrictPhys` with the only difference that damping in normal direction can be considered.",
		((shared_ptr<MatchMaker>,kn,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's normal contact stiffnesses. If this value is not given the elastic properties (i.e. young) of the two colliding materials are used to calculate the stiffness."))
		((shared_ptr<MatchMaker>,kRatio,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's shear contact stiffnesses. If this value is not given the elastic properties (i.e. poisson) of the two colliding materials are used to calculate the stiffness."))
		((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's friction angle. If ``None``, minimum value is used."))
		);
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictViscoMat_FrictViscoPhys);

/** 2d functor creating the interaction law (Law2) based on SphereContactGeometry (ScGeom) and FrictViscoPhys of 2 bodies, returning type FrictViscoPM */
class Law2_ScGeom_FrictViscoPhys_CundallStrackVisco: public LawFunctor{
	public:
		OpenMPAccumulator<Real> plasticDissipation;
		virtual bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I);
		Real elasticEnergy ();
		Real getPlasticDissipation();
		void initPlasticDissipation(Real initVal=0);
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGeom_FrictViscoPhys_CundallStrackVisco,LawFunctor,"Constitutive law for the FrictViscoPM. Corresponds to :yref:`Law2_ScGeom_FrictPhys_CundallStrack` with the only difference that viscous damping in normal direction can be considered.",
		((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
		((bool,sphericalBodies,true,,"If true, compute branch vectors from radii (faster), else use contactPoint-position. Turning this flag true is safe for sphere-sphere contacts and a few other specific cases. It will give wrong values of torques on facets or boxes."))
		((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts. This will trace only plastic energy in this law, see O.trackEnergy for a more complete energies tracing"))
		((int,plastDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for plastic dissipation (with O.trackEnergy)"))
		((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
		,,
		.def("elasticEnergy",&Law2_ScGeom_FrictViscoPhys_CundallStrackVisco::elasticEnergy,"Compute and return the total elastic energy in all \"FrictViscoPhys\" contacts")
		.def("plasticDissipation",&Law2_ScGeom_FrictViscoPhys_CundallStrackVisco::getPlasticDissipation,"Total energy dissipated in plastic slips at all FrictPhys contacts. Computed only if :yref:Law2_ScGeom_FrictViscoPhys_CundallStrackVisco::traceEnergy` is true.")
		.def("initPlasticDissipation",&Law2_ScGeom_FrictViscoPhys_CundallStrackVisco::initPlasticDissipation,"Initialize cummulated plastic dissipation to a value (0 by default).")
		);
		FUNCTOR2D(ScGeom,FrictViscoPhys);
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_ScGeom_FrictViscoPhys_CundallStrackVisco);


