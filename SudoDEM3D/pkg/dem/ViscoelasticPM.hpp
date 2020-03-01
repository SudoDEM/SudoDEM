// 2009 © Sergei Dorofeenko <sega@users.berlios.de>
// This file contains a set of classes for modelling of viscoelastic
// particles.

#pragma once

#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/pkg/dem/DemXDofGeom.hpp>
#include<sudodem/pkg/common/MatchMaker.hpp>

#ifdef SUDODEM_SPH
#include<sudodem/pkg/common/SPHEngine.hpp>
#endif


/* Simple viscoelastic model */

/// Material
/// Note: Shop::getViscoelasticFromSpheresInteraction can get kn,cn,ks,cs from a analytical solution of a pair spheres interaction problem.
class ViscElMat : public FrictMat {
	public:
		virtual ~ViscElMat();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(ViscElMat,FrictMat,"Material for simple viscoelastic model of contact from analytical solution of a pair spheres interaction problem  [Pournin2001]_ .",
		((Real,tc,NaN,,"Contact time"))
		((Real,en,NaN,,"Restitution coefficient in normal direction"))
		((Real,et,NaN,,"Restitution coefficient in tangential direction"))
		((Real,kn,NaN,,"Normal elastic stiffness. Attention, this parameter cannot be set if tc, en or es is defined!"))
		((Real,cn,NaN,,"Normal viscous constant. Attention, this parameter cannot be set if tc, en or es is defined!"))
		((Real,ks,NaN,,"Shear elastic stiffness. Attention, this parameter cannot be set if tc, en or es is defined!"))
		((Real,cs,NaN,,"Shear viscous constant. Attention, this parameter cannot be set if tc, en or es is defined!"))
		((Real,mR,0.0,,"Rolling resistance, see [Zhou1999536]_."))
#ifdef SUDODEM_SPH
		((bool,SPHmode,false,,"True, if SPH-mode is enabled."))
		((Real,mu,-1,, "Viscosity. See Mueller [Mueller2003]_ ."))                                              // [Mueller2003], (14)
		((int,KernFunctionPressure,Spiky,, "Kernel function for pressure calculation (by default - Spiky). The following kernel functions are available: Poly6=1, Spiky=2, Visco=3, Lucy=4, Monaghan=5."))
		((int,KernFunctionVisco,   Visco,, "Kernel function for viscosity calculation (by default - Visco). The following kernel functions are available: Poly6=1, Spiky=2, Visco=3, Lucy=4, Monaghan=5."))
#endif
		((unsigned int,mRtype,1,,"Rolling resistance type, see [Zhou1999536]_. mRtype=1 - equation (3) in [Zhou1999536]_; mRtype=2 - equation (4) in [Zhou1999536]_.")),
		createIndex();
	);
	REGISTER_CLASS_INDEX(ViscElMat,FrictMat);
};
REGISTER_SERIALIZABLE(ViscElMat);

/// Interaction physics
class ViscElPhys : public FrictPhys{
	public:
		virtual ~ViscElPhys();
		Real R;
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(ViscElPhys,FrictPhys,"IPhys created from :yref:`ViscElMat`, for use with :yref:`Law2_ScGeom_ViscElPhys_Basic`.",
		((Real,cn,NaN,,"Normal viscous constant"))
		((Real,cs,NaN,,"Shear viscous constant"))
		((Real,mR,0.0,,"Rolling resistance, see [Zhou1999536]_."))
#ifdef SUDODEM_SPH
		((bool,SPHmode,false,,"True, if SPH-mode is enabled."))
		((Real,h,-1,,    "Core radius. See Mueller [Mueller2003]_ ."))                                            // [Mueller2003], (1)
		((Real,mu,-1,,   "Viscosity. See Mueller [Mueller2003]_ ."))                                              // [Mueller2003], (14)
#endif
		((unsigned int,mRtype,1,,"Rolling resistance type, see [Zhou1999536]_. mRtype=1 - equation (3) in [Zhou1999536]_; mRtype=2 - equation (4) in [Zhou1999536]_")),
		createIndex();
	)
#ifdef SUDODEM_SPH
		KernelFunction kernelFunctionCurrentPressure;
		KernelFunction kernelFunctionCurrentVisco;
#endif
	REGISTER_CLASS_INDEX(ViscElPhys,FrictPhys);
};
REGISTER_SERIALIZABLE(ViscElPhys);

/// Convert material to interaction physics.
// Uses the rule of consecutively connection.
class Ip2_ViscElMat_ViscElMat_ViscElPhys: public IPhysFunctor {
	public :
		virtual void go(const shared_ptr<Material>& b1,
					const shared_ptr<Material>& b2,
					const shared_ptr<Interaction>& interaction);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_ViscElMat_ViscElMat_ViscElPhys,IPhysFunctor,"Convert 2 instances of :yref:`ViscElMat` to :yref:`ViscElPhys` using the rule of consecutive connection.",
 		((shared_ptr<MatchMaker>,tc,,,"Instance of :yref:`MatchMaker` determining contact time"))
		((shared_ptr<MatchMaker>,en,,,"Instance of :yref:`MatchMaker` determining restitution coefficient in normal direction"))
		((shared_ptr<MatchMaker>,et,,,"Instance of :yref:`MatchMaker` determining restitution coefficient in tangential direction")));
	virtual void Calculate_ViscElMat_ViscElMat_ViscElPhys(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction, shared_ptr<ViscElPhys> phys);
	FUNCTOR2D(ViscElMat,ViscElMat);
};
REGISTER_SERIALIZABLE(Ip2_ViscElMat_ViscElMat_ViscElPhys);

/// Constitutive law
/// This class provides linear viscoelastic contact model
class Law2_ScGeom_ViscElPhys_Basic: public LawFunctor {
	public :
		virtual bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
	public :
	FUNCTOR2D(ScGeom,ViscElPhys);
	SUDODEM_CLASS_BASE_DOC(Law2_ScGeom_ViscElPhys_Basic,LawFunctor,"Linear viscoelastic model operating on :yref:`ScGeom` and :yref:`ViscElPhys`. The model is mostly based on the paper for For details see Pournin [Pournin2001]_ .");
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_ScGeom_ViscElPhys_Basic);

Real contactParameterCalculation(const Real& l1,const Real& l2);
bool computeForceTorqueViscEl(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I, Vector3r & force, Vector3r & torque1, Vector3r & torque2);
