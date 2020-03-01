#include<sudodem/core/Clump.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/dem/Integrator.hpp>

#ifdef SUDODEM_OPENMP
  #include<omp.h>
#endif

SUDODEM_PLUGIN((Integrator));



void observer::operator()( const stateVector&  x , Real  t ) const
{
	this->integrator->scene->time=t;

	this->integrator->setCurrentStates(x);

}


//! Integrator's pseudo-ctor (factory), taking nested lists of slave engines (might be moved to real ctor perhaps)


void Integrator::action(){


}

void Integrator::system(const stateVector& currentstates, stateVector& derivatives, Real time)
{

	#ifdef SUDODEM_OPENMP
	//prevent https://bugs.launchpad.net/sudodem/+bug/923929
		ensureSync();
	#endif

	//Calculate orientation

	maxVelocitySq=-1;

	setCurrentStates(currentstates);

	scene->time=time;

	const int size=(int)slaves.size();

	for(int i=0; i<size; i++){
		// run every slave group sequentially
		FOREACH(const shared_ptr<Engine>& e, slaves[i]) {
			e->scene=scene;
			if(!e->dead && e->isActivated()) e->action();
		}
	}
	derivatives=getSceneStateDot();

/*
	std::cout<<std::endl<<"Derivatives are"<<std::endl;
	for(long int k=0;k<derivatives.size();k++)
	std::cout<<std::endl<<derivatives[k]<<std::endl;*/
}

stateVector& Integrator::getSceneStateDot(void){

	try{

		long int numberofscenebodies=scene->bodies->size();

		scene->forces.sync();

		accumstatedotofthescene.resize(2*scene->bodies->size()*7);

		SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){

		const Body::id_t& id=b->getId();

	 	Vector3r force=Vector3r::Zero();

		Vector3r vel_current;

		Vector3r moment=Vector3r::Zero();

		Vector3r angvel_current;

		Quaternionr ori_current;

		Quaternionr angvelquat;

		Quaternionr oridot_current;

		if(!b->isClumpMember()) {

			Real mass=b->state->mass;

			Vector3r inertia=b->state->inertia;

			vel_current=b->state->vel;

			angvel_current=b->state->angVel;

			ori_current=b->state->ori;

			// clumps forces
			if(b->isClump()) {
				b->shape->cast<Clump>().addForceTorqueFromMembers(b->state.get(),scene,force,moment);
				#ifdef SUDODEM_OPENMP
				//it is safe here, since only one thread will read/write
				scene->forces.addTorqueUnsynced(id,moment);
				scene->forces.addForceUnsynced(id,force);
				#else
				scene->forces.addTorque(id,moment);
				scene->forces.addForce(id,force);
				#endif
			}


			force=scene->forces.getForce(id); moment=scene->forces.getTorque(id);

			/*
			 *	Calculation of accelerations
			 *
			*/


			force[0]=force[0]/mass;	force[1]= force[1]/mass;	force[2]= force[2]/mass; //Calculate linear acceleration

			moment[0]=moment[0]/inertia[0];	moment[1]= moment[1]/inertia[1];	moment[2]= moment[2]/inertia[2]; //Calculate angular acceleration

			//Check for fixation
			/*
				This code block needs optimization
			*/
			string str="xyzXYZ"; // Very very very hard coding!!!!! Fixation seems not handled fine by state structure, should be improved.

			for(int i=0; i<3; i++) if(b->state->blockedDOFs_vec_get().find(str[i]) != std::string::npos){ force[i]=0;vel_current[i]=0;}
			for(int i=3; i<6; i++) if(b->state->blockedDOFs_vec_get().find(str[i]) != std::string::npos){ moment[i-3]=0;angvel_current[i-3]=0;}

			angvelquat = Quaternionr(0.0,angvel_current[0],angvel_current[1],angvel_current[2]);

			oridot_current=0.5*angvelquat*ori_current;


			//	if (densityScaling) accel*=state->densityScaling;

		}
		else
		{

			//is clump member
			force=Vector3r::Zero();
			moment=Vector3r::Zero();
			vel_current=Vector3r::Zero();
		        angvel_current=Vector3r::Zero();
			oridot_current=Quaternionr(0,0,0,0);// zero change in quaternion with respect to time
		}


		/*Orientation differantion is straight forward.*/
		accumstatedotofthescene[id*7+0]=vel_current[0]; 	accumstatedotofthescene[(id+numberofscenebodies)*7+0]=force[0];
		accumstatedotofthescene[id*7+1]=vel_current[1];		accumstatedotofthescene[(id+numberofscenebodies)*7+1]=force[1];
		accumstatedotofthescene[id*7+2]=vel_current[2];		accumstatedotofthescene[(id+numberofscenebodies)*7+2]=force[2];
		accumstatedotofthescene[id*7+3]=oridot_current.w();	accumstatedotofthescene[(id+numberofscenebodies)*7+3]=moment[0];
		accumstatedotofthescene[id*7+4]=oridot_current.x();	accumstatedotofthescene[(id+numberofscenebodies)*7+4]=moment[1];
		accumstatedotofthescene[id*7+5]=oridot_current.y();	accumstatedotofthescene[(id+numberofscenebodies)*7+5]=moment[2];
		accumstatedotofthescene[id*7+6]=oridot_current.z();	accumstatedotofthescene[(id+numberofscenebodies)*7+6]=0;


		} SUDODEM_PARALLEL_FOREACH_BODY_END();


	}
	catch(std::exception& e){

		LOG_FATAL("Unhandled exception at Integrator::getSceneStateDot the exception information : "<<typeid(e).name()<<" : "<<e.what());
	}

	return accumstatedotofthescene;



}



stateVector& Integrator::getCurrentStates(void)
{


	try{

		long int numberofscenebodies=scene->bodies->size();

		accumstateofthescene.resize(2*scene->bodies->size()*7);

		SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){

		const Body::id_t& id=b->getId();

         	Vector3r pos_current=b->state->pos;

		Vector3r vel_current=b->state->vel;

		Quaternionr ori=b->state->ori;

		Vector3r angvel=b->state->angVel;

		accumstateofthescene[id*7+0]=pos_current[0]; 	accumstateofthescene[(id+numberofscenebodies)*7+0]=vel_current[0];
		accumstateofthescene[id*7+1]=pos_current[1];	accumstateofthescene[(id+numberofscenebodies)*7+1]=vel_current[1];
		accumstateofthescene[id*7+2]=pos_current[2];	accumstateofthescene[(id+numberofscenebodies)*7+2]=vel_current[2];
		accumstateofthescene[id*7+3]=ori.w();		accumstateofthescene[(id+numberofscenebodies)*7+3]=angvel[0];
		accumstateofthescene[id*7+4]=ori.x();		accumstateofthescene[(id+numberofscenebodies)*7+4]=angvel[1];
		accumstateofthescene[id*7+5]=ori.y();		accumstateofthescene[(id+numberofscenebodies)*7+5]=angvel[2];
		accumstateofthescene[id*7+6]=ori.z();		accumstateofthescene[(id+numberofscenebodies)*7+6]=0;


		} SUDODEM_PARALLEL_FOREACH_BODY_END();


	}
	catch(std::exception& e){

		LOG_FATAL("Unhandled exception at Integrator::getCurrentStates the exception information : "<<typeid(e).name()<<" : "<<e.what());
	}

	return accumstateofthescene;

}

bool Integrator::setCurrentStates(stateVector yscene)
{

	try{

		long int numberofscenebodies=scene->bodies->size();

		//Zero max velocity for each thread
		#ifdef SUDODEM_OPENMP
			FOREACH(Real& thrMaxVSq, threadMaxVelocitySq) { thrMaxVSq=0; }
		#endif

		SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){

		if(b->isClumpMember()) continue;

		const Body::id_t& id=b->getId();

         	Vector3r pos_current;

		pos_current<<yscene[id*7+0],yscene[id*7+1],yscene[id*7+2];

		Vector3r vel_current;

		vel_current<<yscene[(id+numberofscenebodies)*7+0],yscene[(id+numberofscenebodies)*7+1],yscene[(id+numberofscenebodies)*7+2];

		Quaternionr ori=Quaternionr(yscene[id*7+3],yscene[id*7+4],yscene[id*7+5],yscene[id*7+6]);

		Vector3r angvel;

		angvel<<yscene[(id+numberofscenebodies)*7+3],yscene[(id+numberofscenebodies)*7+4],yscene[(id+numberofscenebodies)*7+5];

                b->state->pos=pos_current;

                b->state->vel=vel_current;

		b->state->ori=ori;

		//std::cout<<"Setting orientation to "<<ori<<std::endl;

		b->state->ori.normalize(); //Normalize orientation

		//std::cout<<"Setting angvel to "<<angvel<<std::endl;

                b->state->angVel=angvel;

		#ifdef SUDODEM_OPENMP
			Real& thrMaxVSq=threadMaxVelocitySq[omp_get_thread_num()]; thrMaxVSq=max(thrMaxVSq,b->state->vel.squaredNorm());
		#else
			maxVelocitySq=max(maxVelocitySq,b->state->vel.squaredNorm());// Set maximum velocity of the scene
		#endif

		if(b->isClump()) Clump::moveMembers(b,scene,this);

		} SUDODEM_PARALLEL_FOREACH_BODY_END();

		#ifdef SUDODEM_OPENMP
			FOREACH(const Real& thrMaxVSq, threadMaxVelocitySq) { maxVelocitySq=max(maxVelocitySq,thrMaxVSq); }
		#endif


	}
	catch(std::exception& e){

		LOG_FATAL("Unhandled exception at Integrator::setCurrentStates the exception information : "<<typeid(e).name()<<" : "<<e.what());

		return false;
	}

	return true;
}


#ifdef SUDODEM_OPENMP
void Integrator::ensureSync()
{
	if (syncEnsured) return;
	SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){
// 		if(b->isClump()) continue;
		scene->forces.addForce(b->getId(),Vector3r(0,0,0));
	} SUDODEM_PARALLEL_FOREACH_BODY_END();
	syncEnsured=true;
}
#endif


void Integrator::saveMaximaDisplacement(const shared_ptr<Body>& b){
	if (!b->bound) return;//clumps for instance, have no bounds, hence not saved
	Vector3r disp=b->state->pos-b->bound->refPos;
	Real maxDisp=max(std::abs(disp[0]),max(std::abs(disp[1]),std::abs(disp[2])));
	if (!maxDisp || maxDisp<b->bound->sweepLength) {/*b->bound->isBounding = (updatingDispFactor>0 && (updatingDispFactor*maxDisp)<b->bound->sweepLength);*/
	maxDisp=0.5;//not 0, else it will be seen as "not updated" by the collider, but less than 1 means no colliding
	}
	else {/*b->bound->isBounding = false;*/ maxDisp=2;/*2 is more than 1, enough to trigger collider*/}

	maxVelocitySq=max(maxVelocitySq,maxDisp);
}

void Integrator::slaves_set(const boost::python::list& slaves2){
std::cout<<"Adding slaves";
	int len=boost::python::len(slaves2);
	slaves.clear();
	for(int i=0; i<len; i++){
		boost::python::extract<std::vector<shared_ptr<Engine> > > serialGroup(slaves2[i]);
		if (serialGroup.check()){ slaves.push_back(serialGroup()); continue; }
		boost::python::extract<shared_ptr<Engine> > serialAlone(slaves2[i]);
		if (serialAlone.check()){ vector<shared_ptr<Engine> > aloneWrap; aloneWrap.push_back(serialAlone()); slaves.push_back(aloneWrap); continue; }
		PyErr_SetString(PyExc_TypeError,"Engines that are given to Integrator should be in two cases (a) in an ordered group, (b) alone engines");
		boost::python::throw_error_already_set();
	}
}

boost::python::list Integrator::slaves_get(){
	boost::python::list ret;
	FOREACH(vector<shared_ptr<Engine > >& grp, slaves){
		if(grp.size()==1) ret.append(boost::python::object(grp[0]));
		else ret.append(boost::python::object(grp));
	}
	return ret;
}


