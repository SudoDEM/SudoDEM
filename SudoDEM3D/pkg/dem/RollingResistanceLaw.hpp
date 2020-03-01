
#pragma once

#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<boost/tuple/tuple.hpp>
#include<sudodem/pkg/common/NormShearPhys.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>


// The following code was moved from RolFrictMat.hpp
class RolFrictMat : public FrictMat
{
	public :
		virtual ~RolFrictMat () {};

/// Serialization
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(RolFrictMat,FrictMat,"",
		//((bool,isCohesive,true,,""))
		((Real,Kn,1e7,,"Normal contact stiffness [N/m]."))
		((Real,Ks,1e7,,"Tangential contact stiffness [N/m]."))
		((Real,alphaKr,2.0,,"Dimensionless rolling stiffness."))
		((Real,alphaKtw,2.0,,"Dimensionless twist stiffness."))
		((Real,etaRoll,-1.,,"Dimensionless rolling (aka 'bending') strength. If negative, rolling moment will be elastic."))
		((Real,etaTwist,-1.,,"Dimensionless twisting strength. If negative, twist moment will be elastic."))
		,
		createIndex();
		);
/// Indexable
	REGISTER_CLASS_INDEX(RolFrictMat,FrictMat);
};

REGISTER_SERIALIZABLE(RolFrictMat);

// The following code was moved from RolFrictPhys.hpp
class RolFrictPhys : public FrictPhys
{
	public :
		virtual ~RolFrictPhys() {};
		//void SetBreakingState() {cohesionBroken = true; normalAdhesion = 0; shearAdhesion = 0;};

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(RolFrictPhys,FrictPhys,"",
		((Real,kr,0,,"rotational stiffness [N.m/rad]"))
		((Real,ktw,0,,"twist stiffness [N.m/rad]"))
		((Real,maxRollPl,0.0,,"Coefficient of rolling friction (negative means elastic)."))
		((Real,maxTwistPl,0.0,,"Coefficient of twisting friction (negative means elastic)."))
		((Real,creep_viscosity,-1,,"creep viscosity [Pa.s/m]."))
		// internal attributes
		((Vector3r,moment_twist,Vector3r(0,0,0),(Attr::noSave | Attr::readonly),"Twist moment"))
		((Vector3r,moment_bending,Vector3r(0,0,0),(Attr::noSave | Attr::readonly),"Bending moment"))
		,
		createIndex();
	);
/// Indexable
	REGISTER_CLASS_INDEX(RolFrictPhys,FrictPhys);
};

REGISTER_SERIALIZABLE(RolFrictPhys);

class RollingResistanceLaw: public LawFunctor{
	public:
		OpenMPAccumulator<Real> plasticDissipation;
		Real normElastEnergy();
		Real shearElastEnergy();
		Real bendingElastEnergy();
		Real twistElastEnergy();
		Real totalElastEnergy();
		Real getPlasticDissipation();
		void initPlasticDissipation(Real initVal=0);
	virtual bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I);
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(RollingResistanceLaw,LawFunctor,"Law for linear traction-compression-bending-twisting, with cohesion+friction and Mohr-Coulomb plasticity surface. This law adds adhesion and moments to :yref:`Law2_ScGeom_FrictPhys_CundallStrack`.\n\nThe normal force is (with the convention of positive tensile forces) $F_n=min(k_n*u_n, a_n)$, with $a_n$ the normal adhesion. The shear force is $F_s=k_s*u_s$, the plasticity condition defines the maximum value of the shear force, by default $F_s^{max}=F_n*tan(\\phi)+a_s$, with $\\phi$ the friction angle and $a_s$ the shear adhesion. If :yref:`RolFrictPhys::cohesionDisableFriction` is True, friction is ignored as long as adhesion is active, and the maximum shear force is only $F_s^{max}=a_s$.\n\nIf the maximum tensile or maximum shear force is reached and :yref:`RolFrictPhys::fragile` =True (default), the cohesive link is broken, and $a_n, a_s$ are set back to zero. If a tensile force is present, the contact is lost, else the shear strength is $F_s^{max}=F_n*tan(\\phi)$. If :yref:`RolFrictPhys::fragile` =False (in course of implementation), the behaviour is perfectly plastic, and the shear strength is kept constant.\n\nIf :yref:`RollingResistanceLaw::momentRotationLaw` =True, bending and twisting moments are computed using a linear law with moduli respectively $k_t$ and $k_r$ (the two values are the same currently), so that the moments are : $M_b=k_b*\\Theta_b$ and $M_t=k_t*\\Theta_t$, with $\\Theta_{b,t}$ the relative rotations between interacting bodies. The maximum value of moments can be defined and takes the form of rolling friction. Cohesive -type moment may also be included in the future.\n\nCreep at contact is implemented in this law, as defined in [Hassan2010]_. If activated, there is a viscous behaviour of the shear and twisting components, and the evolution of the elastic parts of shear displacement and relative twist is given by $du_{s,e}/dt=-F_s/\\nu_s$ and $d\\Theta_{t,e}/dt=-M_t/\\nu_t$.",
		((bool,use_rolling_resistance,true,,"If true, use bending/twisting moments (rolling resistance) at all contacts. "))
		((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts. This will trace only plastic energy in this law, see O.trackEnergy for a more complete energies tracing"))
		((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
		((int,shearDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for shear dissipation (with O.trackEnergy)"))
		((int,bendingDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for bending dissipation (with O.trackEnergy)"))
		((int,twistDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for twist dissipation (with O.trackEnergy)"))
		,,
		.def("normElastEnergy",&RollingResistanceLaw::normElastEnergy,"Compute normal elastic energy.")
		.def("shearElastEnergy",&RollingResistanceLaw::shearElastEnergy,"Compute shear elastic energy.")
		.def("bendingElastEnergy",&RollingResistanceLaw::bendingElastEnergy,"Compute bending elastic energy.")
		.def("twistElastEnergy",&RollingResistanceLaw::twistElastEnergy,"Compute twist elastic energy.")
		.def("elasticEnergy",&RollingResistanceLaw::totalElastEnergy,"Compute total elastic energy.")
	);
	FUNCTOR2D(ScGeom,RolFrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(RollingResistanceLaw);


// The following code was moved from Ip2_RolFrictMat_RolFrictMat_RolFrictPhys.hpp
class Ip2_RolFrictMat_RolFrictMat_RolFrictPhys : public IPhysFunctor
{
	public :
		virtual void go(	const shared_ptr<Material>& b1,
					const shared_ptr<Material>& b2,
					const shared_ptr<Interaction>& interaction);
		

		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(Ip2_RolFrictMat_RolFrictMat_RolFrictPhys,IPhysFunctor,
		"Generates cohesive-frictional interactions with moments, used in the contact law :yref:`RollingResistanceLaw`. The normal/shear stiffness and friction definitions are the same as in :yref:`Ip2_FrictMat_FrictMat_FrictPhys`, check the documentation there for details.\n\nAdhesions related to the normal and the shear components are calculated form :yref:`RolFrictMat::normalCohesion` ($C_n$) and :yref:`RolFrictMat::shearlCohesion` ($C_s$). For particles of size $R_1$,$R_2$ the adhesion will be $a_i=C_i min(R_1,R_2)^2$, $i=n\\,s$.\n\nTwist and rolling stiffnesses are proportional to the shear stiffness through dimensionless factors alphaKtw and alphaKr, such that the rotational stiffnesses are defined by $k_s \\alpha_i R_1 R_2$, $i=tw\\,r$",
		,
		
		);
	FUNCTOR2D(RolFrictMat,RolFrictMat);
};

REGISTER_SERIALIZABLE(Ip2_RolFrictMat_RolFrictMat_RolFrictPhys);
