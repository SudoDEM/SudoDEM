
#include<sudodem/pkg/dem/SpheresFactory.hpp>
#include<sudodem/pkg/common/Sphere.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>



SUDODEM_PLUGIN((SpheresFactory)(CircularFactory)(BoxFactory));
CREATE_LOGGER(SpheresFactory);
CREATE_LOGGER(CircularFactory);
CREATE_LOGGER(BoxFactory);

// initialize random number generator with time seed
static boost::minstd_rand randGen(TimingInfo::getNow(/* get the number even if timing is disabled globally */ true));
static boost::variate_generator<boost::minstd_rand&, boost::uniform_real<Real> > randomUnit(randGen, boost::uniform_real<Real>(0,1));

void SpheresFactory::pickRandomPosition(Vector3r&,Real){
	LOG_FATAL("Engine "<<getClassName()<<" calling virtual method SpheresFactory::pickRandomPosition(), but had to call derived class. This could occur if you use SpheresFactory directly instead derived engines. If not, please submit bug report at http://bugs.launchpad.net/sudodem.");
	throw std::logic_error("SpheresFactory::pickRandomPosition() called.");
}

void SpheresFactory::action(){

	if(!collider){
		FOREACH(const shared_ptr<Engine>& e, scene->engines){ collider=SUDODEM_PTR_DYN_CAST<Collider>(e); if(collider) break; }
		if(!collider) throw runtime_error("SpheresFactory: No Collider instance found in engines (needed for collision detection).");
	}

	goalMass+=massFlowRate*scene->dt; // totalMass that we want to attain in the current step

	if ((PSDcum.size()>0) and (!PSDuse)) {			//Defined, that we will use PSD

		if ((PSDcum.size() != PSDsizes.size()) and (exactDiam)) {								//The number of elements in both arrays should be the same
			LOG_ERROR("PSDcum and PSDsizes should have an equal number of elements, if exactDiam=True.");
			throw std::logic_error("PSDcum and PSDsizes should have an equal number of elements, if exactDiam=True.");
		} else if ((PSDcum.size() != (PSDsizes.size()-1)) and (not(exactDiam))) {//The number of elements in PSDsizes should be on 1 more, than in PSDcum
			LOG_ERROR("PSDsizes should have a number of elements on 1 more, than PSDcum, if exactDiam=False");
			throw std::logic_error("PSDsizes should have a number of elements on 1 more, than PSDcum, if exactDiam=False");
		}

		//Check the correctness of inputted data PSDcum
		for (unsigned int i=1; i<PSDcum.size(); i++) {
			if (PSDcum[i]<PSDcum[i-1] or PSDcum[i-1]<=0) {
				LOG_ERROR("PSDcum should have an ascending positive series of numbers (for example: 0.1, 0.3, 0.5, 1.0)");
				throw std::logic_error("PSDcum should have an ascending positive series of numbers (for example: 0.1, 0.3, 0.5, 1.0)");
			}
		}
		//Check the correctness of inputted data PSDsizes
		for (unsigned int i=1; i<PSDsizes.size(); i++) {
			if (PSDsizes[i]<PSDsizes[i-1] or PSDsizes[i-1]<=0) {
				LOG_ERROR("PSDsizes should have an ascending positive series of numbers (for example: 15, 20, 50, 80)");
				throw std::logic_error("PSDsizes should have an ascending positive series of numbers (for example: 15, 20, 50, 80)");
			}
		}

		//Make normalization of PSDcum
		if (PSDcum[PSDcum.size()]!=1.0) {
			Real k;
			k = 1.0/PSDcum[PSDcum.size()-1];
			for (unsigned int i=1; i<PSDcum.size(); i++) {
				PSDcum[i] = PSDcum[i]*k;
			}
		}

		PSDuse = true;

		//Prepare main vectors
		for (unsigned int i=0; i<PSDcum.size(); i++) {
			if (i==0) {
				PSDNeedProc.push_back(PSDcum[i]);
			} else {
				PSDNeedProc.push_back(PSDcum[i]-PSDcum[i-1]);
			}
			PSDCurMean.push_back(0);
			PSDCurProc.push_back(0);
		}
	}

	normal.normalize();

	LOG_TRACE("totalMass="<<totalMass<<", goalMass="<<goalMass);

	if (maxMass>0 && maxParticles>0) {
		LOG_WARN("Both maxMass and maxParticles cannot be > 0; Setting maxMass=-1)");
		maxMass = -1;
	}

	vector< SpherCoord > justCreatedBodies;

	// pick random initial velocity
	Vector3r initVel;
	if (normalVel.norm()>0) {
		normalVel.normalize();
		initVel = normalVel;
	} else {
		initVel = normal;
	}
	initVel*=(vMin+randomUnit()*(vMax-vMin)); // TODO: compute from vMin, vMax, vAngle, normal;

	while(totalMass<goalMass && (maxParticles<0 || numParticles<maxParticles) && (maxMass<0 || totalMass<maxMass)){
		Real r=0.0;

		Real maxdiff=0.0;		//This and next variable are for PSD-distribution
		int maxdiffID=0;

		if (PSDuse) {
			//find in what "bin" we have maximal difference between required number of material and current:
			for (unsigned int k=0; k<PSDcum.size(); k++) {
				if ((maxdiff < (PSDNeedProc[k]-PSDCurProc[k]))) {
					maxdiff = PSDNeedProc[k]-PSDCurProc[k];
					maxdiffID = k;
				}
			}

			if (exactDiam) {							//Choose the exact diameter
				r=PSDsizes[maxdiffID]/2.0;
			} else {											//Choose the diameter from the range
				Real rMinE = PSDsizes[maxdiffID]/2.0;
				Real rMaxE = PSDsizes[maxdiffID+1]/2.0;
				r=rMinE+randomUnit()*(rMaxE-rMinE);
			}
		} else {
			// pick random radius
			r=rMin+randomUnit()*(rMax-rMin);
		}

		LOG_TRACE("Radius "<<r);
		Vector3r c=Vector3r::Zero();
		// until there is no overlap, pick a random position for the new particle
		int attempt;
		for(attempt=0; attempt<maxAttempt; attempt++){
			pickRandomPosition(c,r);
			LOG_TRACE("Center "<<c);
			Bound b; b.min=c-Vector3r(r,r,r); b.max=c+Vector3r(r,r,r);
			vector<Body::id_t> collidingParticles=collider->probeBoundingVolume(b);   //Check, whether newly created sphere collides with existing bodies

			bool collideWithNewBodies = false;
			if (justCreatedBodies.size()>0) {			//Check, whether newly created sphere collides with bodies from this scope
				for (unsigned int ii = 0; ii < justCreatedBodies.size(); ii++) {
					if ((justCreatedBodies.at(ii).c-c).norm() < (justCreatedBodies.at(ii).r+r)) {	//Bodies intersect
						collideWithNewBodies = true;
						break;
					}
				}
			}

			if(collidingParticles.size()==0 and not(collideWithNewBodies)) break;
			#ifdef SUDODEM_DEBUG
				FOREACH(const Body::id_t& id, collidingParticles) LOG_TRACE(scene->iter<<":"<<attempt<<": collision with #" <<id);
			#endif
		}
		if(attempt==maxAttempt) {
			if (silent) {LOG_INFO("Unable to place new sphere after "<<maxAttempt<<" attempts!");}
			else {LOG_WARN("Unable to place new sphere after "<<maxAttempt<<" attempts!");}
			if (stopIfFailed) {
				massFlowRate=0; goalMass=totalMass;
			}
			return;
		}

		// create particle
		int mId=(materialId>=0 ? materialId : scene->materials.size()+materialId);
		if(mId<0 || (size_t) mId>=scene->materials.size()) throw std::invalid_argument(("SpheresFactory: invalid material id "+boost::lexical_cast<string>(materialId)).c_str());
		const shared_ptr<Material>& material=scene->materials[mId];
		shared_ptr<Body> b(new Body);
		shared_ptr<Sphere> sphere(new Sphere);
		shared_ptr<State> state(material->newAssocState());
		sphere->radius=r;
		state->pos=state->refPos=c;
		if (color[0]>=0 and color[1]>=0 and color[2]>=0){
			sphere->color = color;
		}

		state->vel=initVel;
		Real vol=(4/3.)*Mathr::PI*pow(r,3);
		state->mass=vol*material->density;
		state->inertia=(2./5.)*vol*r*r*material->density*Vector3r::Ones();
		state->blockedDOFs_vec_set(blockedDOFs);

		b->shape=sphere;
		b->state=state;
		b->material=material;
		if (mask>0) {b->groupMask=mask;}
		// insert particle in the simulation
		scene->bodies->insert(b);
		ids.push_back(b->getId());
		// increment total mass, volume and numparticles we've spit out
		totalMass+=state->mass;
		totalVolume+= vol;
		numParticles++;

		justCreatedBodies.push_back(SpherCoord(c, r));

		if (PSDuse) {		//Add newly created "material" into the bin
			Real summMaterial = 0.0;
			if (PSDcalculateMass) { PSDCurMean[maxdiffID]=PSDCurMean[maxdiffID]+state->mass; summMaterial = totalMass;}
			else { PSDCurMean[maxdiffID]=PSDCurMean[maxdiffID]+1; summMaterial = numParticles;}

			for (unsigned int k=0; k<PSDcum.size(); k++) {			//Update  relationships in bins
				PSDCurProc[k] = PSDCurMean[k]/summMaterial;
			}
		}

	}
	//std::cout<<"mass flow rate: "<<totalMass<<endl;totalMass =0.0;
};

void CircularFactory::pickRandomPosition(Vector3r& c, Real r){
	const Quaternionr q(Quaternionr().setFromTwoVectors(Vector3r::UnitZ(),normal));
	Real angle=randomUnit()*2*Mathr::PI, rr=randomUnit()*(radius-r); // random polar coordinate inside the circle
	Real l=(randomUnit()-0.5)*length;
	c = center+q*Vector3r(cos(angle)*rr,sin(angle)*rr,0)+normal*l;
}

void BoxFactory::pickRandomPosition(Vector3r& c, Real r){
	const Quaternionr q(Quaternionr().setFromTwoVectors(Vector3r::UnitZ(),normal));
	//c=center+q*Vector3r((randomUnit()-.5)*2*(extents[0]-r),(randomUnit()-.5)*2*(extents[1]-r),(randomUnit()-.5)*2*(extents[2]-r));
	c=center+q*Vector3r((randomUnit()-.5)*2*(extents[0]),(randomUnit()-.5)*2*(extents[1]),(randomUnit()-.5)*2*(extents[2]));
}
