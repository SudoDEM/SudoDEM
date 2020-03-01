#pragma once

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
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/lib/base/openmp-accu.hpp>

/*! Elastic material */
class VolumeFricMat: public Material{
	public:
		 VolumeFricMat(double N, double S, double F){Kn=N; Ks=S; frictionAngle=F;};
		 double GetStrength(){return strength;};
	virtual ~VolumeFricMat(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(VolumeFricMat,Material,"Elastic material with Coulomb friction.",
		((Real,Kn,1e8,,"Normal volumetric 'stiffness' (N/m3)."))
		((Real,Ks,1e5,,"Shear stiffness (N/m)."))
		((Real,frictionAngle,.5,,"Contact friction angle (in radians)."))
		((bool,IsSplitable,0,,"To be splitted ... or not"))
		((double,strength,100,,"Stress at whis polyhedra of volume 4/3*pi [mm] breaks.")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(VolumeFricMat,Material);
};
REGISTER_SERIALIZABLE(VolumeFricMat);

class VolumeFricPhys: public IPhys{
	public:
	virtual ~VolumeFricPhys(){};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(VolumeFricPhys,IPhys,"Simple elastic material with friction for volumetric constitutive laws",
		((Real,kn,0,,"Normal stiffness"))
		((Vector3r,normalForce,Vector3r::Zero(),,"Normal force after previous step (in global coordinates)."))
		((Real,ks,0,,"Shear stiffness"))
		((Vector3r,shearForce,Vector3r::Zero(),,"Shear force after previous step (in global coordinates)."))
		((Real,tangensOfFrictionAngle,NaN,,"tangens of angle of internal friction")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(VolumeFricPhys,IPhys);
};
REGISTER_SERIALIZABLE(VolumeFricPhys);

//***************************************************************************
class Ip2_VolumeFricMat_VolumeFricMat_VolumeFricPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
			const shared_ptr<Material>& b2,
			const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(VolumeFricMat,VolumeFricMat);
	SUDODEM_CLASS_BASE_DOC_ATTRS(Ip2_VolumeFricMat_VolumeFricMat_VolumeFricPhys,IPhysFunctor,"",
	);
};
REGISTER_SERIALIZABLE(Ip2_VolumeFricMat_VolumeFricMat_VolumeFricPhys);


class VolumetricLaw: public LawFunctor{
	OpenMPAccumulator<Real> plasticDissipation;
	virtual bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
	Real elasticEnergy ();
	Real getPlasticDissipation();
	void initPlasticDissipation(Real initVal=0);
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(VolumetricLaw,LawFunctor,"Calculate physical response of 2 :yref:`vector<Polyhedra>` in interaction, based on penetration configuration given by :yref:`PolyhedraGeom`.",
	((Vector3r,shearForce,Vector3r::Zero(),,"Shear force from last step"))
	((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
	((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts. This will trace only plastic energy in this law, see O.trackEnergy for a more complete energies tracing"))
	((int,plastDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for plastic dissipation (with O.trackEnergy)"))
	((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
	,,
	//.def("elasticEnergy",&VolumeFricVolumetricLaw::elasticEnergy,"Compute and return the total elastic energy in all \"FrictPhys\" contacts")
	//.def("plasticDissipation",&VolumeFricVolumetricLaw::getPlasticDissipation,"Total energy dissipated in plastic slips at all FrictPhys contacts. Computed only if :yref:`Law2_ScGeom_FrictPhys_CundallStrack::traceEnergy` is true.")
	//.def("initPlasticDissipation",&VolumeFricVolumetricLaw::initPlasticDissipation,"Initialize cummulated plastic dissipation to a value (0 by default).")
	);
	FUNCTOR2D(ScGeom,VolumeFricPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(VolumetricLaw);


