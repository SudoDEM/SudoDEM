/*************************************************************************
 Copyright (C) 2017 by Sway Zhao		                                 *
*  zhswee@gmail.com      				                            	 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/fem/FEMdriver.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Clump.hpp>
//#include<sudodem/lib/base/Math.hpp>


SUDODEM_PLUGIN((FEMdriver));
CREATE_LOGGER(FEMdriver);

// 1st order numerical damping
void FEMdriver::cundallDamp1st(Vector3r& force, const Vector3r& vel){
	for(int i=0; i<3; i++) force[i]*=1-damping*Mathr::Sign(force[i]*vel[i]);
}
// 2nd order numerical damping
void FEMdriver::cundallDamp2nd(const Real& dt, const Vector3r& vel, Vector3r& accel){
	for(int i=0; i<3; i++) accel[i]*= 1 - damping*Mathr::Sign ( accel[i]*(vel[i] + 0.5*dt*accel[i]) );
}

Vector3r FEMdriver::computeAccel(const Vector3r& force, const Real& mass, int blockedDOFs){
	if(blockedDOFs==0) return (force/mass + gravity);
	Vector3r ret(Vector3r::Zero());
	for(int i=0; i<3; i++) if(!(blockedDOFs & Node::axisDOF(i,false))) ret[i]+=force[i]/mass+gravity[i];
	return ret;
}
Vector3r FEMdriver::computeAngAccel(const Vector3r& torque, const Vector3r& inertia, int blockedDOFs){
	if(blockedDOFs==0) return torque.cwiseQuotient(inertia);
	Vector3r ret(Vector3r::Zero());
	for(int i=0; i<3; i++) if(!(blockedDOFs & Node::axisDOF(i,true))) ret[i]+=torque[i]/inertia[i];
	return ret;
}

#ifdef SUDODEM_OPENMP
void FEMdriver::ensureSync()
{
	if (syncEnsured) return;
	SUDODEM_PARALLEL_FOREACH_NODE_BEGIN(const shared_ptr<Node>& b, scene->nodes){
// 		if(b->isClump()) continue;
		scene->nodeforces.addForce(b->getId(),Vector3r(0,0,0));
	} SUDODEM_PARALLEL_FOREACH_NODE_END();
	syncEnsured=true;
}
#endif


void FEMdriver::applyNodalForces(){
    Real dt = scene->dt;
    scene->forces.sync();
	SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){
    // clump members are handled inside clumps
    //if(b->isClumpMember()) continue;
    if (b->shape->getClassName()!="TriElement") continue;
     const Body::id_t& id=b->getId();
    const shared_ptr<TriElement>& Elm = SUDODEM_PTR_CAST<TriElement> (b->shape);
    b->state->pos = Elm->pos;
    //b->state->ori = Elm->ori;
	//if(!particle->contacts.empty()) distributeForces(particle,sh->cast<Facet>(),bending||applyBary);
    //
    //force and torque at contact
	Vector3r fc=Vector3r::Zero();
	Vector3r mc=Vector3r::Zero();
    fc=scene->forces.getForce(id); mc=scene->forces.getTorque(id);
	Vector3r weights;
    //calculate the weight to distribute the contact force and toqure on each node
	if(false){
		// find barycentric coordinates of the projected contact point, and use those as weights
		// I *guess* this would be the solution for linear interpolation anyway
		//const Vector3r& c=C->geom->node->pos;
		//Vector3r p=c-(c-f.nodes[0]->pos).dot(normal)*normal;
		//weights=CompUtils::triangleBarycentrics(p,f.nodes[0]->pos,f.nodes[1]->pos,f.nodes[2]->pos);
	} else {
		// distribute equally when bending is not considered
		weights=Vector3r::Constant(1/3.);
	}

    //
	Elm->stepUpdate(dt,rotIncr);
	// assemble local stiffness matrix, in case it does not exist yet
	Elm->ensureStiffnessMatrices(young,nu,thickness,bending,bendThickness);
	// compute nodal forces response here
	// ?? CST forces are applied with the - sign, DKT with the + sign; are uXy/phiXy introduced differently?
	Vector6r Fcst=-(Elm->KKcst*Elm->uXy).transpose();


	Vector9r Fdkt;
	if(bending){
		#ifdef MEMBRANE_CONDENSE_DKT
			assert(Elm->KKdkt.size()==54);
			Fdkt=(Elm->KKdkt*Elm->phiXy).transpose();
		#else
			assert(Elm->KKdkt.size()==81);
			Vector9r uDkt_;
			uDkt_<<0,Elm->phiXy.segment<2>(0),0,Elm->phiXy.segment<2>(2),0,Elm->phiXy.segment<2>(4);
			Fdkt=(Elm->KKdkt*uDkt_).transpose();
			#ifdef MEMBRANE_DEBUG_ROT
				Elm->uDkt=uDkt_; // debugging copy, acessible from python
			#endif
		#endif
	} else {
		Fdkt=Vector9r::Zero();
	}
	LOG_TRACE("CST: "<<Fcst.transpose())
	LOG_TRACE("DKT: "<<Fdkt.transpose())
	// surface load, if any
	Real surfLoadForce=0.;
	if(!isnan(Elm->surfLoad) && Elm->surfLoad!=0.){ surfLoadForce=(1/3.)*Elm->getArea()*Elm->surfLoad; }
	// apply nodal forces
	for(int i:{0,1,2}){
		Vector3r Fl=Vector3r(Fcst[2*i],Fcst[2*i+1],Fdkt[3*i]+surfLoadForce);
		Vector3r Tl=Vector3r(Fdkt[3*i+1],Fdkt[3*i+2],0);
        //scene->nodeforces.addForce(Elm->nodes[i]->getId(),Vector3r(1,1,1));//test openMP
		scene->nodeforces.addForce(Elm->nodes[i]->getId(),Elm->ori*Fl+fc*weights[i]);
		scene->nodeforces.addTorque(Elm->nodes[i]->getId(),Elm->ori*Tl+mc*weights[i]);//FIXME:not distributing torques from nodal force (contact)
		LOG_TRACE("  "<<i<<" F: "<<Fl.transpose()<<" \t| "<<Elm->ori*Fl);
		LOG_TRACE("  "<<i<<" T: "<<Tl.transpose()<<" \t| "<<Elm->ori*Tl);
	}
	} SUDODEM_PARALLEL_FOREACH_BODY_END();
}

void FEMdriver::driveNodes(){
    Real dt = scene->dt;
	SUDODEM_PARALLEL_FOREACH_NODE_BEGIN(const shared_ptr<Node>& b, scene->nodes){

			const Node::id_t& id=b->getId();
			Vector3r f=Vector3r::Zero();
			Vector3r m=Vector3r::Zero();

			//in most cases, the initial force on clumps will be zero and next line is not changing f and m, but make sure we don't miss something (e.g. user defined forces on clumps)
			f=scene->nodeforces.getForce(id); m=scene->nodeforces.getTorque(id);
			#ifdef SUDODEM_DEBUG
				if(isnan(f[0])||isnan(f[1])||isnan(f[2])) throw runtime_error(("FEMdriver: NaN force acting on #"+boost::lexical_cast<string>(id)+".").c_str());
				if(isnan(m[0])||isnan(m[1])||isnan(m[2])) throw runtime_error(("FEMdriver: NaN torque acting on #"+boost::lexical_cast<string>(id)+".").c_str());
				if(b->mass<=0 && ((b->blockedDOFs & Node::DOF_XYZ) != Node::DOF_XYZ)) throw runtime_error(("FEMdriver: #"+boost::lexical_cast<string>(id)+" has some linear accelerations enabled, but State::mass is non-positive."));
				if(b->inertia.minCoeff()<=0 && ((b->blockedDOFs & Node::DOF_RXRYRZ) != Node::DOF_RXRYRZ)) throw runtime_error(("FEMdriver: #"+boost::lexical_cast<string>(id)+" has some angular accelerations enabled, but inertia contains non-positive terms."));
			#endif

			// whether to use aspherical rotation integration for this body; for no accelerations, spherical integrator is "exact" (and faster)
			//bool useAspherical=(exactAsphericalRot && b->isAspherical() && state->blockedDOFs!=State::DOF_ALL);

			// for particles not totally blocked, compute accelerations; otherwise, the computations would be useless
			if (b->blockedDOFs!=State::DOF_ALL) {
				// linear acceleration
				Vector3r linAccel=computeAccel(f,b->mass,b->blockedDOFs);
				//if (densityScaling) linAccel*=state->densityScaling;
				if(b->isDamped) cundallDamp2nd(dt,b->vel,linAccel);
				b->vel+=dt*linAccel;
				// angular acceleration
				//cout<<"useAs="<<useAspherical<<endl;
				
				Vector3r angAccel=computeAngAccel(m,b->inertia,b->blockedDOFs);
				//if (densityScaling) angAccel*=state->densityScaling;
				if(b->isDamped) cundallDamp2nd(dt,b->angVel,angAccel);
				b->angVel+=dt*angAccel;
				
			} 
			// update positions from velocities (or torque, for the aspherical integrator)
			//check quiet_system_flag
            //rotate
            Real angle2=b->angVel.squaredNorm();
	        if (angle2!=0 and ( (mask<=0) or ((mask>0) and (Node::byId(id)->maskCompatible(mask))) )) {//If we have an angular velocity, we make a rotation
		        Real angle=sqrt(angle2);
		        Quaternionr q(AngleAxisr(angle*dt,b->angVel/angle));
		        b->ori = q*b->ori;
	        }
	        if(scene->nodeforces.getMoveRotUsed() && scene->nodeforces.getRot(id)!=Vector3r::Zero()
		        and ( (mask<=0) or ((mask>0) and (Node::byId(id)->maskCompatible(mask))) )) {
		        Vector3r r(scene->nodeforces.getRot(id));
		        Real norm=r.norm(); r/=norm;
		        Quaternionr q(AngleAxisr(norm,r));
		        b->ori=q*b->ori;
	        }
	        b->ori.normalize();
			//leapfrogSphericalRotate(state,id,dt);
			if (scene->nodeforces.getMoveRotUsed()) b->pos+=scene->nodeforces.getMove(id);
	        if ( (mask<=0) or ((mask>0) and (Node::byId(id)->maskCompatible(mask))) ) {
		        b->pos+=b->vel*dt;
	        }    
                //leapfrogTranslate(state,id,dt);
	} SUDODEM_PARALLEL_FOREACH_NODE_END();

}

void FEMdriver::action()
{
	#ifdef SUDODEM_OPENMP
	//prevent https://bugs.launchpad.net/sudodem/+bug/923929
	ensureSync();
	#endif
    scene->nodeforces.reset(scene->iter);

    //apply forces nodes
    applyNodalForces();
	scene->nodeforces.sync();
 
    //drive nodes
    driveNodes();
}


// http://www.euclideanspace.com/physics/kinematics/angularvelocity/QuaternionDifferentiation2.pdf
Quaternionr FEMdriver::DotQ(const Vector3r& angVel, const Quaternionr& Q){
	Quaternionr dotQ;
	dotQ.w() = (-Q.x()*angVel[0]-Q.y()*angVel[1]-Q.z()*angVel[2])/2;
	dotQ.x() = ( Q.w()*angVel[0]-Q.z()*angVel[1]+Q.y()*angVel[2])/2;
	dotQ.y() = ( Q.z()*angVel[0]+Q.w()*angVel[1]-Q.x()*angVel[2])/2;
	dotQ.z() = (-Q.y()*angVel[0]+Q.x()*angVel[1]+Q.w()*angVel[2])/2;
	return dotQ;
}
