
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/KinematicEngines.hpp>
#include<sudodem/pkg/dem/Shop.hpp>
#include<sudodem/lib/smoothing/LinearInterpolate.hpp>

SUDODEM_PLUGIN((KinematicEngine)(CombinedKinematicEngine)(TranslationEngine)(HarmonicMotionEngine)(RotationEngine)(HelixEngine)(InterpolatingHelixEngine)(HarmonicRotationEngine)(ServoPIDController)(BicyclePedalEngine));

CREATE_LOGGER(KinematicEngine);

void KinematicEngine::action(){
	if (ids.size()>0) {
		FOREACH(Body::id_t id,ids){
			assert(id<(Body::id_t)scene->bodies->size());
			Body* b=Body::byId(id,scene).get();
			if(b) b->state->vel=b->state->angVel=Vector3r::Zero();
		}
		apply(ids);
	} else {
		LOG_WARN("The list of ids is empty! Can't move any body.");
	}
}

void CombinedKinematicEngine::action(){
	if (ids.size()>0) {
		// reset first
		FOREACH(Body::id_t id,ids){
			assert(id<(Body::id_t)scene->bodies->size());
			Body* b=Body::byId(id,scene).get();
			if(b) b->state->vel=b->state->angVel=Vector3r::Zero();
		}
		// apply one engine after another
		FOREACH(const shared_ptr<KinematicEngine>& e, comb){
			if (e->dead) continue;
			e->scene=scene; e->apply(ids);
		}
	} else {
		LOG_WARN("The list of ids is empty! Can't move any body.");
	}
}

const shared_ptr<CombinedKinematicEngine> CombinedKinematicEngine::fromTwo(const shared_ptr<KinematicEngine>& first, const shared_ptr<KinematicEngine>& second){
	shared_ptr<CombinedKinematicEngine> ret(new CombinedKinematicEngine);
	ret->ids=first->ids; ret->comb.push_back(first); ret->comb.push_back(second);
	return ret;
}


void TranslationEngine::apply(const vector<Body::id_t>& ids){
	if (ids.size()>0) {
		#ifdef SUDODEM_OPENMP
		const long size=ids.size();
		#pragma omp parallel for schedule(static)
		for(long i=0; i<size; i++){
			const Body::id_t& id=ids[i];
		#else
		FOREACH(Body::id_t id,ids){
		#endif
			assert(id<(Body::id_t)scene->bodies->size());
			Body* b=Body::byId(id,scene).get();
			if(!b) continue;
			b->state->vel+=velocity*translationAxis;
		}
	} else {
		LOG_WARN("The list of ids is empty! Can't move any body.");
	}
}

void HarmonicMotionEngine::apply(const vector<Body::id_t>& ids){
	if (ids.size()>0) {
		Vector3r w = f*2.0*Mathr::PI; 										 								//Angular frequency
		Vector3r velocity = (((w*scene->time + fi).array().sin())*(-1.0));
    velocity = velocity.cwiseProduct(A);
    velocity = velocity.cwiseProduct(w);
		FOREACH(Body::id_t id,ids){
			assert(id<(Body::id_t)scene->bodies->size());
			Body* b=Body::byId(id,scene).get();
			if(!b) continue;
			b->state->vel+=velocity;
		}
	} else {
		LOG_WARN("The list of ids is empty! Can't move any body.");
	}
}


void InterpolatingHelixEngine::apply(const vector<Body::id_t>& ids){
	Real virtTime=wrap ? Shop::periodicWrap(scene->time,*times.begin(),*times.rbegin()) : scene->time;
	angularVelocity=linearInterpolate<Real,Real>(virtTime,times,angularVelocities,_pos);
	linearVelocity=angularVelocity*slope;
	HelixEngine::apply(ids);
}

void HelixEngine::apply(const vector<Body::id_t>& ids){
	if (ids.size()>0) {
		const Real& dt=scene->dt;
		angleTurned+=angularVelocity*dt;
		shared_ptr<BodyContainer> bodies = scene->bodies;
		FOREACH(Body::id_t id,ids){
			assert(id<(Body::id_t)bodies->size());
			Body* b=Body::byId(id,scene).get();
			if(!b) continue;
			b->state->vel+=linearVelocity*rotationAxis;
		}
		rotateAroundZero=true;
		RotationEngine::apply(ids);
	} else {
		LOG_WARN("The list of ids is empty! Can't move any body.");
	}
}

void RotationEngine::apply(const vector<Body::id_t>& ids){
	if (ids.size()>0) {
		#ifdef SUDODEM_OPENMP
		const long size=ids.size();
		#pragma omp parallel for schedule(static)
		for(long i=0; i<size; i++){
			const Body::id_t& id=ids[i];
		#else
		FOREACH(Body::id_t id,ids){
		#endif
			assert(id<(Body::id_t)scene->bodies->size());
			Body* b=Body::byId(id,scene).get();
			if(!b) continue;
			b->state->angVel+=rotationAxis*angularVelocity;
			if(rotateAroundZero){
				const Vector3r l=b->state->pos-zeroPoint;
				Quaternionr q(AngleAxisr(angularVelocity*scene->dt,rotationAxis));
				Vector3r newPos=q*l+zeroPoint;
				b->state->vel+=Vector3r(newPos-b->state->pos)/scene->dt;
			}
		}
	} else {
		LOG_WARN("The list of ids is empty! Can't move any body.");
	}
}

void HarmonicRotationEngine::apply(const vector<Body::id_t>& ids){
	const Real& time=scene->time;
	Real w = f*2.0*Mathr::PI; 			//Angular frequency
	angularVelocity = -1.0*A*w*sin(w*time + fi);
	RotationEngine::apply(ids);
}

void ServoPIDController::apply(const vector<Body::id_t>& ids){

  if (iterPrevStart < 0 or ((scene->iter-iterPrevStart)>=iterPeriod)) {

    Vector3r tmpForce = Vector3r::Zero();

    if (ids.size()>0) {
      FOREACH(Body::id_t id,ids){
        assert(id<(Body::id_t)scene->bodies->size());
        tmpForce += scene->forces.getForce(id);
      }
    } else {
      LOG_WARN("The list of ids is empty!");
    }

    axis.normalize();
    tmpForce = tmpForce.cwiseProduct(axis);   // Take into account given axis
    errorCur = tmpForce.norm() - target;      // Find error

    const Real pTerm = errorCur*kP;               // Calculate proportional term
    iTerm += errorCur*kI;                         // Calculate integral term
    const Real dTerm = (errorCur-errorPrev)*kD;   // Calculate derivative term

    errorPrev = errorCur;                         // Save the current value of the error

    curVel = (pTerm + iTerm + dTerm);             // Calculate current velocity

    if (std::abs(curVel) > std::abs(maxVelocity)) {
      curVel*=std::abs(maxVelocity)/std::abs(curVel);
    }

    iterPrevStart = scene->iter;
    current = tmpForce;
  }

  translationAxis = axis;
  velocity = curVel;

  TranslationEngine::apply(ids);
}


void BicyclePedalEngine::apply(const vector<Body::id_t>& ids){
	if (ids.size()>0) {
		Quaternionr qRotateZVec(Quaternionr().setFromTwoVectors(Vector3r(0,0,1), rotationAxis));

		Vector3r newPos = Vector3r(cos(fi + angularVelocity*scene->dt)*radius, sin(fi + angularVelocity*scene->dt)*radius, 0.0);
		Vector3r oldPos = Vector3r(cos(fi)*radius, sin(fi)*radius, 0.0);

		Vector3r newVel = (oldPos - newPos)/scene->dt;

		fi += angularVelocity*scene->dt;
		newVel =  qRotateZVec*newVel;

		#ifdef SUDODEM_OPENMP
		const long size=ids.size();
		#pragma omp parallel for schedule(static)
		for(long i=0; i<size; i++){
			const Body::id_t& id=ids[i];
		#else
		FOREACH(Body::id_t id,ids){
		#endif
			assert(id<(Body::id_t)scene->bodies->size());
			Body* b=Body::byId(id,scene).get();
			if(!b) continue;
			b->state->vel+=newVel;
		}
	} else {
		LOG_WARN("The list of ids is empty! Can't move any body.");
	}
}
