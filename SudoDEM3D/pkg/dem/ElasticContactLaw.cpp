/*************************************************************************
*  Copyright (C) 2005 by Bruno Chareyre   bruno.chareyre@hmg.inpg.fr     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"ElasticContactLaw.hpp"
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>
#include<sudodem/pkg/dem/DemXDofGeom.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>

SUDODEM_PLUGIN((Law2_ScGeom_FrictPhys_CundallStrack)(Law2_ScGeom_ViscoFrictPhys_CundallStrack)(ElasticContactLaw));

#if 1
Real Law2_ScGeom_FrictPhys_CundallStrack::getPlasticDissipation() {return (Real) plasticDissipation;}
void Law2_ScGeom_FrictPhys_CundallStrack::initPlasticDissipation(Real initVal) {plasticDissipation.reset(); plasticDissipation+=initVal;}
Real Law2_ScGeom_FrictPhys_CundallStrack::elasticEnergy()
{
	Real energy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		FrictPhys* phys = dynamic_cast<FrictPhys*>(I->phys.get());
		if(phys) {
			energy += 0.5*(phys->normalForce.squaredNorm()/phys->kn + phys->shearForce.squaredNorm()/phys->ks);}
	}
	return energy;
}
#endif

void ElasticContactLaw::action()
{
	if(!functor) functor=shared_ptr<Law2_ScGeom_FrictPhys_CundallStrack>(new Law2_ScGeom_FrictPhys_CundallStrack);
	functor->neverErase=neverErase;
	functor->scene=scene;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		#ifdef SUDODEM_DEBUG
			// these checks would be redundant in the functor (LawDispatcher does that already)
			if(!dynamic_cast<ScGeom*>(I->geom.get()) || !dynamic_cast<FrictPhys*>(I->phys.get())) continue;
		#endif
			functor->go(I->geom, I->phys, I.get());
	}
}

CREATE_LOGGER(Law2_ScGeom_FrictPhys_CundallStrack);
bool Law2_ScGeom_FrictPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
	int id1 = contact->getId1(), id2 = contact->getId2();

	ScGeom*    geom= static_cast<ScGeom*>(ig.get());
	FrictPhys* phys = static_cast<FrictPhys*>(ip.get());
	if(geom->penetrationDepth <0){
		if (neverErase) {
			phys->shearForce = Vector3r::Zero();
			phys->normalForce = Vector3r::Zero();}
		else return false;
	}
	Real& un=geom->penetrationDepth;
	phys->normalForce=phys->kn*std::max(un,(Real) 0)*geom->normal;

	Vector3r& shearForce = geom->rotate(phys->shearForce);
	const Vector3r& shearDisp = geom->shearIncrement();
	shearForce -= phys->ks*shearDisp;
	Real maxFs = phys->normalForce.squaredNorm()*std::pow(phys->tangensOfFrictionAngle,2);

	if (!scene->trackEnergy  && !traceEnergy){//Update force but don't compute energy terms (see below))
		// PFC3d SlipModel, is using friction angle. CoulombCriterion
		if( shearForce.squaredNorm() > maxFs ){
			Real ratio = sqrt(maxFs) / shearForce.norm();
			shearForce *= ratio;}
	} else {
		//almost the same with additional Vector3r instatinated for energy tracing,
		//duplicated block to make sure there is no cost for the instanciation of the vector when traceEnergy==false
		if(shearForce.squaredNorm() > maxFs){
			Real ratio = sqrt(maxFs) / shearForce.norm();
			Vector3r trialForce=shearForce;//store prev force for definition of plastic slip
			//define the plastic work input and increment the total plastic energy dissipated
			shearForce *= ratio;
			Real dissip=((1/phys->ks)*(trialForce-shearForce))/*plastic disp*/ .dot(shearForce)/*active force*/;
			if (traceEnergy) plasticDissipation += dissip;
			else if(dissip>0) scene->energy->add(dissip,"plastDissip",plastDissipIx,/*reset*/false);
		}
		// compute elastic energy as well
		scene->energy->add(0.5*(phys->normalForce.squaredNorm()/phys->kn+phys->shearForce.squaredNorm()/phys->ks),"elastPotential",elastPotentialIx,/*reset at every timestep*/true);
	}
	if (!scene->isPeriodic && !sphericalBodies) {
		State* de1 = Body::byId(id1,scene)->state.get();
		State* de2 = Body::byId(id2,scene)->state.get();
		applyForceAtContactPoint(-phys->normalForce-shearForce, geom->contactPoint, id1, de1->se3.position, id2, de2->se3.position);}
	else {//we need to use correct branches in the periodic case, the following apply for spheres only
		Vector3r force = -phys->normalForce-shearForce;
		scene->forces.addForce(id1,force);
		scene->forces.addForce(id2,-force);
		scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(force));
		scene->forces.addTorque(id2,(geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(force));
	}
	return true;
}

bool Law2_ScGeom_ViscoFrictPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
	ScGeom*    geom= static_cast<ScGeom*>(ig.get());
	ViscoFrictPhys* phys = static_cast<ViscoFrictPhys*>(ip.get());
	if (shearCreep) {
			const Real& dt=scene->dt;
			geom->rotate(phys->creepedShear);
			phys->creepedShear+= creepStiffness*phys->ks*(phys->shearForce-phys->creepedShear)*dt/viscosity;
			phys->shearForce -= phys->ks*((phys->shearForce-phys->creepedShear)*dt/viscosity);}
	return Law2_ScGeom_FrictPhys_CundallStrack::go(ig,ip,contact);
}
