#include "RollingResistanceLaw.hpp"
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>

SUDODEM_PLUGIN((RollingResistanceLaw)(RolFrictMat)(RolFrictPhys)(Ip2_RolFrictMat_RolFrictMat_RolFrictPhys));
CREATE_LOGGER(RollingResistanceLaw);

Real RollingResistanceLaw::getPlasticDissipation() {return (Real) plasticDissipation;}
void RollingResistanceLaw::initPlasticDissipation(Real initVal) {plasticDissipation.reset(); plasticDissipation+=initVal;}

Real RollingResistanceLaw::normElastEnergy()
{
	Real normEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		RolFrictPhys* phys = SUDODEM_CAST<RolFrictPhys*>(I->phys.get());
		if (phys) {
			normEnergy += 0.5*(phys->normalForce.squaredNorm()/phys->kn);
		}
	}
	return normEnergy;
}
Real RollingResistanceLaw::shearElastEnergy()
{
	Real shearEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		RolFrictPhys* phys = SUDODEM_CAST<RolFrictPhys*>(I->phys.get());
		if (phys) {
			shearEnergy += 0.5*(phys->shearForce.squaredNorm()/phys->ks);
		}
	}
	return shearEnergy;
}

Real RollingResistanceLaw::bendingElastEnergy()
{
	Real bendingEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		RolFrictPhys* phys = SUDODEM_CAST<RolFrictPhys*>(I->phys.get());
		if (phys) {
			bendingEnergy += 0.5*(phys->moment_bending.squaredNorm()/phys->kr);
		}
	}
	return bendingEnergy;
}

Real RollingResistanceLaw::twistElastEnergy()
{
	Real twistEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		RolFrictPhys* phys = SUDODEM_CAST<RolFrictPhys*>(I->phys.get());
		if (phys) {
			twistEnergy += 0.5*(phys->moment_twist.squaredNorm()/phys->ktw);
		}
	}
	return twistEnergy;
}

Real RollingResistanceLaw::totalElastEnergy()
{
	Real totalEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		RolFrictPhys* phys = SUDODEM_CAST<RolFrictPhys*>(I->phys.get());
		if (phys) {
			totalEnergy += 0.5*(phys->normalForce.squaredNorm()/phys->kn);
			totalEnergy += 0.5*(phys->shearForce.squaredNorm()/phys->ks);
			totalEnergy += 0.5*(phys->moment_bending.squaredNorm()/phys->kr);
			totalEnergy += 0.5*(phys->moment_twist.squaredNorm()/phys->ktw);
		}
	}
	return totalEnergy;
}


bool RollingResistanceLaw::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	const Real& dt = scene->dt;
	const int &id1 = contact->getId1();
	const int &id2 = contact->getId2();
	ScGeom* geom  = SUDODEM_CAST<ScGeom*> (ig.get());
	RolFrictPhys* phys = SUDODEM_CAST<RolFrictPhys*> (ip.get());
	Vector3r& shearForce  = phys->shearForce;

	if (contact->isFresh(scene)) shearForce   = Vector3r::Zero();
	Real un     = geom->penetrationDepth;
	//normal contact force
	Real Fn    = phys->kn*std::max(un,(Real) 0);

	phys->normalForce = Fn*geom->normal;
	State* de1 = Body::byId(id1,scene)->state.get();
	State* de2 = Body::byId(id2,scene)->state.get();
	//tangential contact force
	shearForce = geom->rotate(phys->shearForce);
	
	const Vector3r& shearDisp = geom->shearIncrement();
	shearForce -= phys->ks*shearDisp;
	Real maxFs = phys->normalForce.squaredNorm()*std::pow(phys->tangensOfFrictionAngle,2);
	if (!scene->trackEnergy  &&!traceEnergy){//Update force but don't compute energy terms (see below))
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
			else if(dissip>0) scene->energy->add(dissip,"shearDissip",shearDissipIx,/*reset*/false);
		}
		// compute elastic energy as well
		scene->energy->add(0.5*(phys->normalForce.squaredNorm()/phys->kn+phys->shearForce.squaredNorm()/phys->ks),"elastPotential",elastPotentialIx,/*reset at every timestep*/true);
	}

	
	//Apply the force
	applyForceAtContactPoint(-phys->normalForce-shearForce, geom->contactPoint, id1, de1->se3.position, id2, de2->se3.position + (scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist): Vector3r::Zero()));

	/// rolling resistance (including twisting resistance)
	if (use_rolling_resistance) {
	 	// Use incremental formulation to compute moment_twis and moment_bending
		Vector3r relAngVel = geom->getRelAngVel(de1,de2,dt);
		
		//relative angular velocity in the twisting direction
		Vector3r relAngVelTwist = geom->normal.dot(relAngVel)*geom->normal;
		//relative angular velicity in the rolling direction, i.e., bending direction.
		Vector3r relAngVelBend = relAngVel - relAngVelTwist;
		//relative rotational angles due to rolling and twisting
		Vector3r relRotBend = relAngVelBend*dt; //rolling
		Vector3r relRotTwist = relAngVelTwist*dt; // twisting	
		//bending moment and torsional moment
		Vector3r& momentBend = phys->moment_bending;
		Vector3r& momentTwist = phys->moment_twist;
		if (contact->isFresh(scene)){//
			momentBend   = Vector3r::Zero();
			momentTwist   = Vector3r::Zero();
		}
		momentBend = geom->rotate(momentBend); // rotate moment vector (updated)
		momentBend = momentBend-phys->kr*relRotBend;
		momentTwist = geom->rotate(momentTwist); // rotate moment vector (updated)
		momentTwist = momentTwist-phys->ktw*relRotTwist; // FIXME: sign?
		
		/// Plasticity ///
		// limit rolling moment to the plastic value, if required
		if (phys->maxRollPl>=0.){ // do we want to apply plasticity?
			Real RollMax = phys->maxRollPl*phys->normalForce.norm();
			Real scalarRoll = phys->moment_bending.norm();
			if (scalarRoll>RollMax){ // fix maximum rolling moment
				Real ratio = RollMax/scalarRoll;
				phys->moment_bending *= ratio;
				if (scene->trackEnergy){
					Real bendingdissip=((1/phys->kr)*(scalarRoll-RollMax)*RollMax)/*active force*/;
					if(bendingdissip>0) scene->energy->add(bendingdissip,"bendingDissip",bendingDissipIx,/*reset*/false);}
			}
		}
		// limit twisting moment to the plastic value, if required
		if (phys->maxTwistPl>=0.){ // do we want to apply plasticity?
			Real TwistMax = phys->maxTwistPl*phys->normalForce.norm();
			Real scalarTwist= phys->moment_twist.norm();
			if (scalarTwist>TwistMax){ // fix maximum rolling moment
				Real ratio = TwistMax/scalarTwist;
				phys->moment_twist *= ratio;
				if (scene->trackEnergy){
					Real twistdissip=((1/phys->ktw)*(scalarTwist-TwistMax)*TwistMax)/*active force*/;
					if(twistdissip>0) scene->energy->add(twistdissip,"twistDissip",twistDissipIx,/*reset*/false);}
			}	
		}
		// Apply moments now
		Vector3r moment = phys->moment_twist + phys->moment_bending;
		scene->forces.addTorque(id1,-moment);
		scene->forces.addTorque(id2, moment);			
	}
	/// Moment law END       ///
	
	return true;
}


void Ip2_RolFrictMat_RolFrictMat_RolFrictPhys::go(const shared_ptr<Material>& b1    // RolFrictMat
                                        , const shared_ptr<Material>& b2 // RolFrictMat
                                        , const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	RolFrictMat* mat1 = static_cast<RolFrictMat*>(b1.get());
	RolFrictMat* mat2 = static_cast<RolFrictMat*>(b2.get());
	ScGeom* geom = SUDODEM_CAST<ScGeom*>(interaction->geom.get());

	interaction->phys = shared_ptr<RolFrictPhys>(new RolFrictPhys());
	RolFrictPhys* contactPhysics = SUDODEM_CAST<RolFrictPhys*>(interaction->phys.get());

	
    const Body::id_t id= interaction->id1;

   // if (id>5){
        Real frictionAngle = std::min(mat1->frictionAngle,mat2->frictionAngle);
        contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle);
        contactPhysics->kn = 2.0*mat1->Kn*mat2->Kn/(mat1->Kn+mat2->Kn);
	    contactPhysics->ks = 2.0*mat1->Ks*mat2->Ks/(mat1->Ks+mat2->Ks);

		Real Da 	= geom->radius1;
		Real Db 	= geom->radius2;
		// harmonic average of alphas parameters
		Real AlphaKr = 2.0*mat1->alphaKr*mat2->alphaKr/(mat1->alphaKr+mat2->alphaKr);
		Real AlphaKtw;
		if (mat1->alphaKtw && mat2->alphaKtw) AlphaKtw = 2.0*mat1->alphaKtw*mat2->alphaKtw/(mat1->alphaKtw+mat2->alphaKtw);
		else AlphaKtw=0;
		/*
		contactPhysics->kr = Da*Db*contactPhysics->kn*AlphaKr;
		contactPhysics->ktw = Da*Db*contactPhysics->ks*AlphaKtw;
		contactPhysics->maxRollPl = min(mat1->etaRoll*Da,mat2->etaRoll*Db);
		contactPhysics->maxTwistPl = contactPhysics->tangensOfFrictionAngle*min(mat1->etaTwist*Da,mat2->etaTwist*Db);
		*/
		//using a model similar to Jiang et al. (2015)
		contactPhysics->kr = 0.25*Da*Db*contactPhysics->kn*AlphaKr*AlphaKr;
		contactPhysics->ktw = 0.5*Da*Db*contactPhysics->ks*AlphaKtw*AlphaKtw;
		contactPhysics->maxRollPl = 0.25*AlphaKr*min(mat1->etaRoll*Da,mat2->etaRoll*Db);
		contactPhysics->maxTwistPl = 0.65*std::tan(frictionAngle)*AlphaKr*min(Da,Db);
    //}else{//walls, and the wall number is less than 6 by default.
    //    contactPhysics->tangensOfFrictionAngle = std::tan(mat1->frictionAngle);
    //    contactPhysics->kn = mat1->Kn;
	//    contactPhysics->ks = mat1->Ks;
   // }

	//contactPhysics->betan = std::max(mat1->betan,mat2->betan);//
	//contactPhysics->betas = std::max(mat1->betas,mat2->betas);
};

