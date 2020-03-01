/*************************************************************************
 Copyright (C) 2018 by Shiwei Zhao	                         *
*  zhswee@gmail.com      					 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/dem/NewtonIntegrator.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Clump.hpp>
//#include<sudodem/lib/base/Math.hpp>


SUDODEM_PLUGIN((NewtonIntegrator));
CREATE_LOGGER(NewtonIntegrator);

// 1st order numerical damping
void NewtonIntegrator::cundallDamp1st(Vector2r& force, const Vector2r& vel){
	for(int i=0; i<2; i++) force[i]*=1-damping*Mathr::Sign(force[i]*vel[i]);
}
// 2nd order numerical damping
void NewtonIntegrator::cundallDamp2nd(const Real& dt, const Vector2r& vel, Vector2r& accel){
	for(int i=0; i<2; i++) accel[i]*= 1 - damping*Mathr::Sign ( accel[i]*(vel[i] + 0.5*dt*accel[i]) );
}

// 1st order numerical damping
void NewtonIntegrator::cundallDamp1st(Real& m, const Real& angvel){
	m*=1-damping*Mathr::Sign(m*angvel);
}
// 2nd order numerical damping
void NewtonIntegrator::cundallDamp2nd(const Real& dt, const Real& angvel, Real& accel){
	accel*= 1 - damping*Mathr::Sign ( accel*(angvel + 0.5*dt*accel) );//damping for rotation should not be too large, otherwise the torque will become opposite periodically
}

Vector2r NewtonIntegrator::computeAccel(const Vector2r& force, const Real& mass, int blockedDOFs){
	if(blockedDOFs==0) return (force/mass + gravity);
	Vector2r ret(Vector2r::Zero());
	for(int i=0; i<2; i++) if(!(blockedDOFs & State::axisDOF(i,false))) ret[i]+=force[i]/mass+gravity[i];
	return ret;
}
Real NewtonIntegrator::computeAngAccel(const Real& torque, const Real& inertia, int blockedDOFs){
	if(blockedDOFs & State::DOF_RZ){
    return 0.0;
  }else{
		return torque/inertia;
	}
}

void NewtonIntegrator::updateEnergy(const shared_ptr<Body>& b, const State* state, const Vector2r& fluctVel, const Vector2r& f, const Real& m){
	assert(b->isStandalone() || b->isClump());
	// always positive dissipation, by-component: |F_i|*|v_i|*damping*dt (|T_i|*|ω_i|*damping*dt for rotations)
	/*if(damping!=0. && state->isDamped){
		scene->energy->add(fluctVel.cwiseAbs().dot(f.cwiseAbs())*damping*scene->dt,"nonviscDamp",nonviscDampIx,false);
		// when the aspherical integrator is used, torque is damped instead of ang acceleration; this code is only approximate
		scene->energy->add(state->angVel.cwiseAbs().dot(m.cwiseAbs())*damping*scene->dt,"nonviscDamp",nonviscDampIx,false);
	}
	// kinetic energy
	Real Etrans=.5*state->mass*fluctVel.squaredNorm();
	Real Erot;
	// rotational terms
	if(b->isAspherical()){
		Matrix3r mI; mI<<state->inertia[0],0,0, 0,state->inertia[1],0, 0,0,state->inertia[2];
		Matrix3r T(state->ori);
		Erot=.5*b->state->angVel.transpose().dot((T.transpose()*mI*T)*b->state->angVel);
	} else { Erot=0.5*state->angVel.dot(state->inertia.cwiseProduct(state->angVel)); }
	if(!kinSplit) scene->energy->add(Etrans+Erot,"kinetic",kinEnergyIx,true);
	else{ scene->energy->add(Etrans,"kinTrans",kinEnergyTransIx,true); scene->energy->add(Erot,"kinRot",kinEnergyRotIx,true); }
	// gravitational work (work done by gravity is "negative", since the energy appears in the system from outside)
	scene->energy->add(-gravity.dot(b->state->vel)*b->state->mass*scene->dt,"gravWork",fieldWorkIx,false);

*/
}

void NewtonIntegrator::saveMaximaVelocity(const Body::id_t& id, State* state){
	#ifdef SUDODEM_OPENMP
		Real& thrMaxVSq=threadMaxVelocitySq[omp_get_thread_num()]; thrMaxVSq=max(thrMaxVSq,state->vel.squaredNorm());
	#else
		maxVelocitySq=max(maxVelocitySq,state->vel.squaredNorm());
	#endif
}


void NewtonIntegrator::saveMaximaDisplacement(const shared_ptr<Body>& b){
	if (!b->bound) return;//clumps for instance, have no bounds, hence not saved
	Vector2r disp=b->state->pos-b->bound->refPos;
	Real maxDisp=max(std::abs(disp[0]),max(std::abs(disp[1]),std::abs(disp[2])));
	if (!maxDisp || maxDisp<b->bound->sweepLength) {/*b->bound->isBounding = (updatingDispFactor>0 && (updatingDispFactor*maxDisp)<b->bound->sweepLength);*/
	maxDisp=0.5;//not 0, else it will be seen as "not updated" by the collider, but less than 1 means no colliding
	}
	else {/*b->bound->isBounding = false;*/ maxDisp=2;/*2 is more than 1, enough to trigger collider*/}
	#ifdef SUDODEM_OPENMP
		Real& thrMaxVSq=threadMaxVelocitySq[omp_get_thread_num()]; thrMaxVSq=max(thrMaxVSq,maxDisp);
	#else
		maxVelocitySq=max(maxVelocitySq,maxDisp);
	#endif
}

#ifdef SUDODEM_OPENMP
void NewtonIntegrator::ensureSync()
{
	if (syncEnsured) return;
	SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){
// 		if(b->isClump()) continue;
		scene->forces.addForce(b->getId(),Vector2r(0,0));
	} SUDODEM_PARALLEL_FOREACH_BODY_END();
	syncEnsured=true;
}
#endif

void NewtonIntegrator::action()
{
	#ifdef SUDODEM_OPENMP
	//prevent https://bugs.launchpad.net/sudodem/+bug/923929
	ensureSync();
	#endif

	scene->forces.sync();
	bodySelected=(scene->selectedBody>=0);
	if(warnNoForceReset && scene->forces.lastReset<scene->iter) LOG_WARN("O.forces last reset in step "<<scene->forces.lastReset<<", while the current step is "<<scene->iter<<". Did you forget to include ForceResetter in O.engines?");
	const Real& dt=scene->dt;
	//Take care of user's request to change velGrad. Safe to change it here after the interaction loop.
	if (scene->cell->velGradChanged || scene->cell->nextVelGrad!=Matrix2r::Zero()) {
		scene->cell->velGrad=scene->cell->nextVelGrad;
		scene->cell->velGradChanged=0; scene->cell->nextVelGrad=Matrix2r::Zero();}
	homoDeform=scene->cell->homoDeform;
	dVelGrad=scene->cell->velGrad-prevVelGrad;
	//Matrix2r R=.5*(dVelGrad-dVelGrad.transpose());//dSpin= -R(0,1);
	//Real dSpin = -(dVelGrad(0,1) - dVelGrad(1,0));

	//cout<<"dSpin="<<dSpin<<endl;
	// account for motion of the periodic boundary, if we remember its last position
	// its velocity will count as max velocity of bodies
	// otherwise the collider might not run if only the cell were changing without any particle motion
	// FIXME: will not work for pure shear transformation, which does not change Cell::getSize()
	if(scene->isPeriodic && ((prevCellSize!=scene->cell->getSize())) && /* initial value */!isnan(prevCellSize[0]) ){ cellChanged=true; maxVelocitySq=(prevCellSize-scene->cell->getSize()).squaredNorm()/pow(dt,2); }
	else { maxVelocitySq=0; cellChanged=false; }

	#ifdef SUDODEM_BODY_CALLBACK
		// setup callbacks
		vector<BodyCallback::FuncPtr> callbackPtrs;
		FOREACH(const shared_ptr<BodyCallback>& cb, callbacks){
			cerr<<"<cb="<<cb.get()<<", setting cb->scene="<<scene<<">";
			cb->scene=scene;
			callbackPtrs.push_back(cb->stepInit());
		}
		assert(callbackPtrs.size()==callbacks.size());
		size_t callbacksSize=callbacks.size();
	#endif

	const bool trackEnergy(scene->trackEnergy);
	const bool isPeriodic(scene->isPeriodic);

	#ifdef SUDODEM_OPENMP
		FOREACH(Real& thrMaxVSq, threadMaxVelocitySq) { thrMaxVSq=0; }
	#endif
	SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){
			// clump members are handled inside clumps
            if (b->shape->getClassName()=="TriElement") continue;
			if(b->isClumpMember()) continue;
			State* state=b->state.get(); const Body::id_t& id=b->getId();
			Vector2r f=Vector2r::Zero();
			Real m=0.0;

			// clumps forces
			if(b->isClump()) {
				b->shape->cast<Clump>().addForceTorqueFromMembers(state,scene,f,m);
				#ifdef SUDODEM_OPENMP
				//it is safe here, since only one thread is adding forces/torques
				scene->forces.addTorqueUnsynced(id,m);
				scene->forces.addForceUnsynced(id,f);
				#else
				scene->forces.addTorque(id,m);
				scene->forces.addForce(id,f);
				#endif
			}
			//in most cases, the initial force on clumps will be zero and next line is not changing f and m, but make sure we don't miss something (e.g. user defined forces on clumps)
			f=scene->forces.getForce(id); m=scene->forces.getTorque(id);
			#ifdef SUDODEM_DEBUG
				if(isnan(f[0])||isnan(f[1]) throw runtime_error(("NewtonIntegrator: NaN force acting on #"+boost::lexical_cast<string>(id)+".").c_str());
				if(isnan(m)) throw runtime_error(("NewtonIntegrator: NaN torque acting on #"+boost::lexical_cast<string>(id)+".").c_str());
				if(state->mass<=0 && ((state->blockedDOFs & State::DOF_XYZ) != State::DOF_XYZ)) throw runtime_error(("NewtonIntegrator: #"+boost::lexical_cast<string>(id)+" has some linear accelerations enabled, but State::mass is non-positive."));
				if(state->inertia.minCoeff()<=0 && ((state->blockedDOFs & State::DOF_RXRYRZ) != State::DOF_RXRYRZ)) throw runtime_error(("NewtonIntegrator: #"+boost::lexical_cast<string>(id)+" has some angular accelerations enabled, but State::inertia contains non-positive terms."));
			#endif

			// fluctuation velocity does not contain meanfield velocity in periodic boundaries
			// in aperiodic boundaries, it is equal to absolute velocity
			Vector2r fluctVel=isPeriodic?scene->cell->bodyFluctuationVel(b->state->pos,b->state->vel,prevVelGrad):state->vel;

			// numerical damping & kinetic energy
			if(trackEnergy) updateEnergy(b,state,fluctVel,f,m);

			// whether to use aspherical rotation integration for this body; for no accelerations, spherical integrator is "exact" (and faster)
			bool useAspherical=(exactAsphericalRot && b->isAspherical() && state->blockedDOFs!=State::DOF_ALL);

			// for particles not totally blocked, compute accelerations; otherwise, the computations would be useless
			if (state->blockedDOFs!=State::DOF_ALL) {
				// linear acceleration
				Vector2r linAccel=computeAccel(f,state->mass,state->blockedDOFs);
				if (densityScaling) linAccel*=state->densityScaling;
				if(state->isDamped) cundallDamp2nd(dt,fluctVel,linAccel);
				//This is the convective term, appearing in the time derivation of Cundall/Thornton expression (dx/dt=velGrad*pos -> d²x/dt²=dvelGrad/dt*pos+velGrad*vel), negligible in many cases but not for high speed large deformations (gaz or turbulent flow).
				if (isPeriodic && homoDeform) linAccel+=prevVelGrad*state->vel;
				//finally update velocity
				state->vel+=dt*linAccel;
				// angular acceleration
				//cout<<"useAs="<<useAspherical<<endl;
				if(!useAspherical){ // uses angular velocity
					//cout<<"m="<<m<<endl;
					Real angAccel=computeAngAccel(m,state->inertia,state->blockedDOFs);
					if (densityScaling) angAccel*=state->densityScaling;
					if(state->isDamped) cundallDamp2nd(dt,state->angVel,angAccel);
					state->angVel+=dt*angAccel;
					//cout<<"angAccel="<<angAccel<<"anglVel="<<state->angVel<<endl;
				} else { // uses torque
					if(state->blockedDOFs & State::DOF_RZ) {m=0;state->angVel=0.0;} // block DOFs here
					else{
						if(state->isDamped){
							//not use torque now for 2d
							Real angAccel=computeAngAccel(m,state->inertia,state->blockedDOFs);
							if (densityScaling) angAccel*=state->densityScaling;
							if(state->isDamped) cundallDamp2nd(dt,state->angVel,angAccel);
							state->angVel+=dt*angAccel;
						  //cundallDamp1st(m,state->angVel);
						}
					}
				}
			// reflect macro-deformation even for non-dynamic bodies
			} else if (isPeriodic && homoDeform) state->vel+=dt*prevVelGrad*state->vel;
			//if (isPeriodic && homoDeform) state->angVel+= dSpin;
			// update positions from velocities (or torque, for the aspherical integrator)
			//check quiet_system_flag
            if (quiet_system_flag){//the flag is set to true from the other implementation
                state->vel = Vector2r::Zero();
                state->angVel = 0.0;
            }
            else{//execute the normal script
								if (isSuperellipse){
			            Superellipse* A = static_cast<Superellipse*>(b->shape.get());
			            if (A->isSphere){
			                    leapfrogSphericalRotate(state,id,dt);
			            }else{  //cout<<"superellipse rotation"<<endl;
			                    leapfrogSuperellipseRotate(A,state,id,dt);
			            }
			          }else{leapfrogSphericalRotate(state,id,dt);}
                leapfrogTranslate(state,id,dt);
            }
			saveMaximaDisplacement(b);
			// move individual members of the clump, save maxima velocity (for collider stride)
			if(b->isClump()) Clump::moveMembers(b,scene,this);

			#ifdef SUDODEM_BODY_CALLBACK
				// process callbacks
				for(size_t i=0; i<callbacksSize; i++){
					cerr<<"<"<<b->id<<",cb="<<callbacks[i]<<",scene="<<callbacks[i]->scene<<">"; // <<",force="<<callbacks[i]->scene->forces.getForce(b->id)<<">";
					if(callbackPtrs[i]!=NULL) (*(callbackPtrs[i]))(callbacks[i].get(),b.get());
				}
			#endif
	} SUDODEM_PARALLEL_FOREACH_BODY_END();
	if (quiet_system_flag){cout<<"quiet_system!"<<endl;quiet_system_flag = false;}
	#ifdef SUDODEM_OPENMP
		FOREACH(const Real& thrMaxVSq, threadMaxVelocitySq) { maxVelocitySq=max(maxVelocitySq,thrMaxVSq); }
	#endif
	if(scene->isPeriodic) { prevCellSize=scene->cell->getSize(); prevVelGrad=scene->cell->prevVelGrad=scene->cell->velGrad; }
}

void NewtonIntegrator::leapfrogTranslate(State* state, const Body::id_t& id, const Real& dt){
	if (scene->forces.getMoveRotUsed()) state->pos+=scene->forces.getMove(id);
	// update velocity reflecting changes in the macroscopic velocity field, making the problem homothetic.
	//NOTE : if the velocity is updated before moving the body, it means the current velGrad (i.e. before integration in cell->integrateAndUpdate) will be effective for the current time-step. Is it correct? If not, this velocity update can be moved just after "state->pos += state->vel*dt", meaning the current velocity impulse will be applied at next iteration, after the contact law. (All this assuming the ordering is resetForces->integrateAndUpdate->contactLaw->PeriCompressor->NewtonsLaw. Any other might fool us.)
	//NOTE : dVel defined without wraping the coordinates means bodies out of the (0,0,0) period can move realy fast. It has to be compensated properly in the definition of relative velocities (see Ig2 functors and contact laws).
		//Reflect mean-field (periodic cell) acceleration in the velocity
	if(scene->isPeriodic && homoDeform) {Vector2r dVel=dVelGrad*state->pos; state->vel+=dVel;}

	if ( (mask<=0) or ((mask>0) and (Body::byId(id)->maskCompatible(mask))) ) {
		state->pos+=state->vel*dt;
	}
}

void NewtonIntegrator::leapfrogSphericalRotate(State* state, const Body::id_t& id, const Real& dt )
{
	Real angle=state->angVel;
	if (angle!=0 and ( (mask<=0) or ((mask>0) and (Body::byId(id)->maskCompatible(mask))) )) {//If we have an angular velocity, we make a rotation
		//Real angle=sqrt(angle2);
		//Quaternionr q(AngleAxisr(angle*dt,state->angVel/angle));
    Rotationr rot(angle*dt);
		//if(id==399){
		//	cout<<"rot angle="<<rot.angle()<<endl;
		//	cout<<"state angle="<<state->ori.angle()<<endl;
		//}
		state->ori = rot*state->ori;
		//state->ori.angle() = state->ori.smallestAngle();//the rotation angle in [-pi,pi]
		//if(id==399){
		//	cout<<"state angle2="<<state->ori.angle()<<endl;
		//}
	}
	//state->ori.normalize();
}

void NewtonIntegrator::leapfrogSuperellipseRotate(Superellipse* shape,State* state, const Body::id_t& id, const Real& dt){

    //Matrix2r A = shape->rot_mat2local;
		Real angle=state->angVel;
		if (angle!=0 and ( (mask<=0) or ((mask>0) and (Body::byId(id)->maskCompatible(mask))) )) {//If we have an angular velocity, we make a rotation
			//Real angle=sqrt(angle2);
			//Quaternionr q(AngleAxisr(angle*dt,state->angVel/angle));
	    Rotationr rot(angle*dt);
			//if(id==399){
			//	cout<<"rot angle="<<rot.angle()<<endl;
			//	cout<<"state angle="<<state->ori.angle()<<endl;
			//}
			state->ori = rot*state->ori;
			//cout<<"here in NewtonEngine! orientation is"<<state->ori.angle()<<endl;
      shape->rot_mat2local = state->ori.inverse().toRotationMatrix();//to particle's system
      shape->rot_mat2global = state->ori.toRotationMatrix(); //to global system
	  }
}


bool NewtonIntegrator::get_densityScaling() {
	//FOREACH(const shared_ptr<Engine> e, Omega::instance().getScene()->engines) {
		//GlobalStiffnessTimeStepper* ts=dynamic_cast<GlobalStiffnessTimeStepper*>(e.get());
		//if (ts && densityScaling != ts->densityScaling) LOG_WARN("density scaling is not active in the timeStepper, it will have no effect unless a scaling is specified manually for some bodies");}
	LOG_WARN("GlobalStiffnessTimeStepper not present in O.engines, density scaling will have no effect unless a scaling is specified manually for some bodies");
	return 1;//densityScaling;
}

void NewtonIntegrator::set_densityScaling(bool dsc) {
	/*FOREACH(const shared_ptr<Engine> e, Omega::instance().getScene()->engines) {
		GlobalStiffnessTimeStepper* ts=dynamic_cast<GlobalStiffnessTimeStepper*>(e.get());
		if (ts) {
			ts->densityScaling=dsc;
			densityScaling=dsc;
			LOG_WARN("GlobalStiffnessTimeStepper found in O.engines and adjusted to match this setting. Revert in the the timestepper if you don't want the scaling adjusted automatically.");
			return;
		}
	} LOG_WARN("GlobalStiffnessTimeStepper not found in O.engines. Density scaling will have no effect unless a scaling is specified manually for some bodies");

*/
densityScaling=dsc;//not used
}
