// 2009 © Sergei Dorofeenko <sega@users.berlios.de>
#include"ViscoelasticPM.hpp"
#include<sudodem/core/State.hpp>
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/Sphere.hpp>

#ifdef SUDODEM_SPH
#include<sudodem/pkg/common/SPHEngine.hpp>
#endif

using std::isfinite;
SUDODEM_PLUGIN((ViscElMat)(ViscElPhys)(Ip2_ViscElMat_ViscElMat_ViscElPhys)(Law2_ScGeom_ViscElPhys_Basic));

/* ViscElMat */
ViscElMat::~ViscElMat(){}

/* ViscElPhys */
ViscElPhys::~ViscElPhys(){}

/* Ip2_ViscElMat_ViscElMat_ViscElPhys */
void Ip2_ViscElMat_ViscElMat_ViscElPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) {
	// no updates of an existing contact
	if(interaction->phys) return;
	shared_ptr<ViscElPhys> phys (new ViscElPhys());
	Calculate_ViscElMat_ViscElMat_ViscElPhys(b1, b2, interaction, phys);
	interaction->phys = phys;
}

/* Law2_ScGeom_ViscElPhys_Basic */
bool Law2_ScGeom_ViscElPhys_Basic::go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) {
	Vector3r force = Vector3r::Zero();
	Vector3r torque1 = Vector3r::Zero();
	Vector3r torque2 = Vector3r::Zero();
	if (computeForceTorqueViscEl(_geom, _phys, I, force, torque1, torque2) and (I->isActive)) {
		const int id1 = I->getId1();
		const int id2 = I->getId2();

		addForce (id1,-force,scene);
		addForce (id2, force,scene);
		addTorque(id1, torque1,scene);
		addTorque(id2, torque2,scene);
		return true;
	} else return false;
}

bool computeForceTorqueViscEl(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I, Vector3r & force, Vector3r & torque1, Vector3r & torque2) {
	ViscElPhys& phys=*static_cast<ViscElPhys*>(_phys.get());
	const ScGeom& geom=*static_cast<ScGeom*>(_geom.get());
	Scene* scene=Omega::instance().getScene().get();

#ifdef SUDODEM_SPH
//=======================================================================================================
	if (phys.SPHmode) {
		if (computeForceSPH(_geom, _phys, I, force)) {
			return true;
		} else {
			return false;
		}
	}
//=======================================================================================================
#endif

	const int id1 = I->getId1();
	const int id2 = I->getId2();

	if (geom.penetrationDepth<0) {
		return false;
	} else {
		const BodyContainer& bodies = *scene->bodies;

		const State& de1 = *static_cast<State*>(bodies[id1]->state.get());
		const State& de2 = *static_cast<State*>(bodies[id2]->state.get());

		Vector3r& shearForce = phys.shearForce;
		if (I->isFresh(scene)) shearForce=Vector3r(0,0,0);
		const Real& dt = scene->dt;
		shearForce = geom.rotate(shearForce);

		// Handle periodicity.
		const Vector3r shift2 = scene->isPeriodic ? scene->cell->intrShiftPos(I->cellDist): Vector3r::Zero();
		const Vector3r shiftVel = scene->isPeriodic ? scene->cell->intrShiftVel(I->cellDist): Vector3r::Zero();

		const Vector3r c1x = (geom.contactPoint - de1.pos);
		const Vector3r c2x = (geom.contactPoint - de2.pos - shift2);

		const Vector3r relativeVelocity = (de1.vel+de1.angVel.cross(c1x)) - (de2.vel+de2.angVel.cross(c2x)) + shiftVel;
		const Real normalVelocity	= geom.normal.dot(relativeVelocity);
		const Vector3r shearVelocity	= relativeVelocity-normalVelocity*geom.normal;

		// As Chiara Modenese suggest, we store the elastic part
		// and then add the viscous part if we pass the Mohr-Coulomb criterion.
		// See http://www.mail-archive.com/sudodem-users@lists.launchpad.net/msg01391.html
		shearForce += phys.ks*dt*shearVelocity; // the elastic shear force have a history, but
		Vector3r shearForceVisc = Vector3r::Zero(); // the viscous shear damping haven't a history because it is a function of the instant velocity


		// Prevent appearing of attraction forces due to a viscous component
		// [Radjai2011], page 3, equation [1.7]
		// [Schwager2007]
		const Real normForceReal = phys.kn * geom.penetrationDepth + phys.cn * normalVelocity;
		if (normForceReal < 0) {
			phys.normalForce = Vector3r::Zero();
		} else {
			phys.normalForce = normForceReal * geom.normal;
		}

		Vector3r momentResistance = Vector3r::Zero();
		if (phys.mR>0.0) {
			const Vector3r relAngVel  = de1.angVel - de2.angVel;
			relAngVel.normalized();

			if (phys.mRtype == 1) {
				momentResistance = -phys.mR*phys.normalForce.norm()*relAngVel;																														// [Zhou1999536], equation (3)
			} else if (phys.mRtype == 2) {
				momentResistance = -phys.mR*(c1x.cross(de1.angVel) - c2x.cross(de2.angVel)).norm()*phys.normalForce.norm()*relAngVel;			// [Zhou1999536], equation (4)
			}
		}

		const Real maxFs = phys.normalForce.squaredNorm() * std::pow(phys.tangensOfFrictionAngle,2);
		if( shearForce.squaredNorm() > maxFs )
		{
			// Then Mohr-Coulomb is violated (so, we slip),
			// we have the max value of the shear force, so
			// we consider only friction damping.
			const Real ratio = sqrt(maxFs) / shearForce.norm();
			shearForce *= ratio;
		}
		else
		{
			// Then no slip occurs we consider friction damping + viscous damping.
			shearForceVisc = phys.cs*shearVelocity;
		}
		force = phys.normalForce + shearForce + shearForceVisc;
		torque1 = -c1x.cross(force)+momentResistance;
		torque2 =  c2x.cross(force)-momentResistance;
		return true;
	}
}

void Ip2_ViscElMat_ViscElMat_ViscElPhys::Calculate_ViscElMat_ViscElMat_ViscElPhys(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction, shared_ptr<ViscElPhys> phys) {
	ViscElMat* mat1 = static_cast<ViscElMat*>(b1.get());
	ViscElMat* mat2 = static_cast<ViscElMat*>(b2.get());
	Real mass1 = 1.0;
	Real mass2 = 1.0;

	if ((isfinite(mat1->kn) and  not (isfinite(mat2->kn))) or
			(isfinite(mat2->kn) and  not (isfinite(mat1->kn))) or
			(isfinite(mat1->ks) and  not (isfinite(mat2->ks))) or
			(isfinite(mat2->ks) and  not (isfinite(mat1->ks))) or
			(isfinite(mat1->cn) and  not (isfinite(mat2->cn))) or
			(isfinite(mat2->cn) and  not (isfinite(mat1->cn))) or
			(isfinite(mat1->cs) and  not (isfinite(mat2->cs))) or
			(isfinite(mat2->cs) and  not (isfinite(mat1->cs))) or
			(isfinite(mat1->tc) and  not (isfinite(mat2->tc))) or
			(isfinite(mat2->tc) and  not (isfinite(mat1->tc))) or
			(isfinite(mat1->en) and  not (isfinite(mat2->en))) or
			(isfinite(mat2->en) and  not (isfinite(mat1->en))) or
			(isfinite(mat1->et) and  not (isfinite(mat2->et))) or
			(isfinite(mat2->et) and  not (isfinite(mat1->et)))) {
				throw runtime_error("Both materials should have the same defined set of variables e.g. tc, ks etc.!");
			}

	mass1 = Body::byId(interaction->getId1())->state->mass;
	mass2 = Body::byId(interaction->getId2())->state->mass;
	if (mass1 == 0.0 and mass2 > 0.0) {
		mass1 = mass2;
	} else if (mass2 == 0.0 and mass1 > 0.0) {
		mass2 = mass1;
	}

	// See [Pournin2001, just below equation (19)]
	const Real massR = mass1*mass2/(mass1+mass2);

	GenericSpheresContact* sphCont=SUDODEM_CAST<GenericSpheresContact*>(interaction->geom.get());
	Real R1=sphCont->refR1>0?sphCont->refR1:sphCont->refR2;
	Real R2=sphCont->refR2>0?sphCont->refR2:sphCont->refR1;

	Real kn1 = 0.0; Real kn2 = 0.0;
	Real cn1 = 0.0; Real cn2 = 0.0;
	Real ks1 = 0.0; Real ks2 = 0.0;
	Real cs1 = 0.0; Real cs2 = 0.0;

	if (((isfinite(mat1->tc)) and (isfinite(mat1->en)) and (isfinite(mat1->et)))  or ((tc) and (en) and (et))) {
		//Set parameters according to [Pournin2001]

		const Real Tc = (tc) ? (*tc)(mat1->id,mat2->id) : (mat1->tc+mat2->tc)/2.0;
		const Real En = (en) ? (*en)(mat1->id,mat2->id) : (mat1->en+mat2->en)/2.0;
		const Real Et = (et) ? (*et)(mat1->id,mat2->id) : (mat1->et+mat2->et)/2.0;

    // Factor 2 at the end of each expression is necessary, because we calculate
    // individual kn1, kn2, ks1, ks2 etc., because kn1 = 2*kn, ks1 = 2*ks
    // http://www.mail-archive.com/sudodem-users@lists.launchpad.net/msg08778.html
    kn1 = kn2 = 1/Tc/Tc * ( Mathr::PI*Mathr::PI + pow(log(En),2) )*massR*2;
    cn1 = cn2 = -2.0 /Tc * log(En)*massR*2;
    ks1 = ks2 = 2.0/7.0 /Tc/Tc * ( Mathr::PI*Mathr::PI + pow(log(Et),2) )*massR*2;
    cs1 = cs2 = -4.0/7.0 /Tc * log(Et)*massR*2;
    //           ^^^
    // It seems to be an error in [Pournin2001] (22) Eq.4, missing factor 2
    // Thanks to Dominik Boemer for pointing this out
    // http://www.mail-archive.com/sudodem-users@lists.launchpad.net/msg08741.html

		if (std::abs(cn1) <= Mathr::ZERO_TOLERANCE ) cn1=0;
		if (std::abs(cn2) <= Mathr::ZERO_TOLERANCE ) cn2=0;
		if (std::abs(cs1) <= Mathr::ZERO_TOLERANCE ) cs1=0;
		if (std::abs(cs2) <= Mathr::ZERO_TOLERANCE ) cs2=0;
	} else if ((isfinite(mat1->kn)) and (isfinite(mat1->ks)) and (isfinite(mat1->cn)) and (isfinite(mat1->cs))) {
		//Set parameters explicitly
		kn1 = mat1->kn;
		kn2 = mat2->kn;
		ks1 = mat1->ks;
		ks2 = mat2->ks;
		cn1 = mat1->cn;
		cn2 = mat2->cn;
		cs1 = mat1->cs;
		cs2 = mat2->cs;
	} else {
		//Set parameters on the base of young modulus
		kn1 = 2*mat1->young*R1;
		kn2 = 2*mat2->young*R2;
		ks1 = kn1*mat1->poisson;
		ks2 = kn2*mat2->poisson;
		if ((isfinite(mat1->cn)) and (isfinite(mat1->cs))) {
			cn1 = mat1->cn;
			cn2 = mat2->cn;
			cs1 = mat1->cs;
			cs2 = mat2->cs;
		}
	}

	const Real mR1 = mat1->mR;      const Real mR2 = mat2->mR;
	const int mRtype1 = mat1->mRtype; const int mRtype2 = mat2->mRtype;


	phys->kn = contactParameterCalculation(kn1,kn2);
	phys->ks = contactParameterCalculation(ks1,ks2);
	phys->cn = contactParameterCalculation(cn1,cn2);
	phys->cs = contactParameterCalculation(cs1,cs2);

 	if ((mR1>0) or (mR2>0)) {
		phys->mR = 2.0/( ((mR1>0)?1/mR1:0) + ((mR2>0)?1/mR2:0) );
	} else {
		phys->mR = 0;
	}

	phys->tangensOfFrictionAngle = std::tan(std::min(mat1->frictionAngle, mat2->frictionAngle));
	phys->shearForce = Vector3r(0,0,0);

	if ((mRtype1 != mRtype2) or (mRtype1>2) or (mRtype2>2) or (mRtype1<1) or (mRtype2<1) ) {
		throw runtime_error("mRtype should be equal for both materials and have the values 1 or 2!");
	} else {
		phys->mRtype = mRtype1;
	}
#ifdef SUDODEM_SPH
	if (mat1->SPHmode and mat2->SPHmode)  {
		phys->SPHmode=true;
		phys->mu=(mat1->mu+mat2->mu)/2.0;
	}

	phys->kernelFunctionCurrentPressure = returnKernelFunction (mat1->KernFunctionPressure, mat2->KernFunctionPressure, Grad);
	phys->kernelFunctionCurrentVisco    = returnKernelFunction (mat1->KernFunctionVisco, mat2->KernFunctionVisco, Lapl);
#endif
}

/* Contact parameter calculation function */
Real contactParameterCalculation(const Real& l1, const Real& l2){
  // If one of paramaters > 0. we DO NOT return 0
  Real a = (l1?1/l1:0) + (l2?1/l2:0);
  if (a) return 1/a;
  else return 0;
}

