#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/State.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>

#include"VolumeFric.hpp"

SUDODEM_PLUGIN(/* self-contained in hpp: */  (VolumeFricPhys) (VolumeFricMat) (Ip2_VolumeFricMat_VolumeFricMat_VolumeFricPhys) (VolumetricLaw)
	/* some code in cpp (this file): */
	);
/* Material law, physics */

void Ip2_VolumeFricMat_VolumeFricMat_VolumeFricPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<VolumeFricMat>& mat1 = SUDODEM_PTR_CAST<VolumeFricMat>(b1);
	const shared_ptr<VolumeFricMat>& mat2 = SUDODEM_PTR_CAST<VolumeFricMat>(b2);
	interaction->phys = shared_ptr<VolumeFricPhys>(new VolumeFricPhys());
	const shared_ptr<VolumeFricPhys>& contactPhysics = SUDODEM_PTR_CAST<VolumeFricPhys>(interaction->phys);
	Real Kna 	= mat1->Kn;
	Real Knb 	= mat2->Kn;
	Real Ksa 	= mat1->Ks;
	Real Ksb 	= mat2->Ks;
	Real frictionAngle = std::min(mat1->frictionAngle,mat2->frictionAngle);
        contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle);
	contactPhysics->kn = Kna*Knb/(Kna+Knb);
	contactPhysics->ks = Ksa*Ksb/(Ksa+Ksb);
};


//**************************************************************************************
// Apply forces on polyhedrons in collision based on geometric configuration
bool VolumetricLaw::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
	int id1 = contact->getId1(), id2 = contact->getId2();

	ScGeom*    geom= static_cast<ScGeom*>(ig.get());
	//FrictPhys* phys = static_cast<FrictPhys*>(ip.get());
	VolumeFricPhys* phys = dynamic_cast<VolumeFricPhys*>(contact->phys.get());
	if(geom->penetrationDepth <0){
		if (neverErase) {
			phys->shearForce = Vector3r::Zero();
			phys->normalForce = Vector3r::Zero();}
		else return false;
	}
	//Real& un=geom->penetrationDepth;
	//calculate normal force in terms of volumetric stiffness
	//calculate the volume of overlap
	//sphere and sphere

	Real ra = geom->radius1;
	Real rb = geom->radius2;
	Real dab = ra + rb - geom->penetrationDepth;
	Real ha = ra - 0.5*dab - (std::pow(ra,2)-std::pow(rb,2))/dab*0.5;
	Real hb = rb - 0.5*dab + (std::pow(ra,2)-std::pow(rb,2))/dab*0.5;
	Real Volume = Mathr::PI/3.*(std::pow(ha,2)*(3.*ra - ha) + std::pow(hb,2)*(3.*rb - hb));
	phys->normalForce=phys->kn*std::max(Volume,(Real) 0)*geom->normal;

	Vector3r& shearForce = geom->rotate(phys->shearForce);
	const Vector3r& shearDisp = geom->shearIncrement();
	shearForce -= phys->ks*shearDisp;
	Real maxFs = phys->normalForce.squaredNorm()*std::pow(phys->tangensOfFrictionAngle,2);
	// PFC3d SlipModel, is using friction angle. CoulombCriterion
	if( shearForce.squaredNorm() > maxFs ){
		Real ratio = sqrt(maxFs) / shearForce.norm();
		shearForce *= ratio;}


	//we need to use correct branches in the periodic case, the following apply for spheres only
	Vector3r force = -phys->normalForce-shearForce;
	scene->forces.addForce(id1,force);
	scene->forces.addForce(id2,-force);
	scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(force));
	scene->forces.addTorque(id2,(geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(force));

	return true;

}
