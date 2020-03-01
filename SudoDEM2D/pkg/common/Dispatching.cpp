#include<sudodem/pkg/common/Dispatching.hpp>

SUDODEM_PLUGIN((BoundFunctor)(IGeomFunctor)(IPhysFunctor)(LawFunctor)(BoundDispatcher)(IGeomDispatcher)(IPhysDispatcher)(LawDispatcher));
BoundFunctor::~BoundFunctor(){};
IGeomFunctor::~IGeomFunctor(){};
IPhysFunctor::~IPhysFunctor(){};
LawFunctor::~LawFunctor(){};


/********************************************************************
                      BoundDispatcher
*********************************************************************/

CREATE_LOGGER(BoundDispatcher);
void BoundDispatcher::action()
{
	updateScenePtr();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	const long numBodies=(long)bodies->size();
	#pragma omp parallel for num_threads(ompThreads>0 ? min(ompThreads,omp_get_max_threads()) : omp_get_max_threads())
	for(int id=0; id<numBodies; id++){
		if(!bodies->exists(id)) continue; // don't delete this check  - Janek
		const shared_ptr<Body>& b=(*bodies)[id];
		processBody(b);
	}
// 	With -j4, this update takes more time that the dispatching in itslef, and it is quite useless: commented out
// 	scene->updateBound();
}

void BoundDispatcher::processBody(const shared_ptr<Body>& b)
{
// 	const shared_ptr<Body>& b=(*bodies)[id];
		shared_ptr<Shape>& shape=b->shape;
		if(!b->isBounded() || !shape) return;
		if(b->bound) {
			Real& sweepLength = b->bound->sweepLength;
			if (targetInterv>=0) {
				Vector2r disp = b->state->pos-b->bound->refPos;
				Real dist = max(std::abs(disp[0]),std::abs(disp[1]));
				if (dist){
					Real newLength = dist*targetInterv/(scene->iter-b->bound->lastUpdateIter);
					newLength = max(0.9*sweepLength,newLength);//don't decrease size too fast to prevent time consuming oscillations
					sweepLength=max(minSweepDistFactor*sweepDist,min(newLength,sweepDist));}
				else sweepLength=0;
			} else sweepLength=sweepDist;
		}
		#ifdef BV_FUNCTOR_CACHE
		if(!shape->boundFunctor){ shape->boundFunctor=this->getFunctor1D(shape); if(!shape->boundFunctor) return; }
		shape->boundFunctor->go(shape,b->bound,b->state->se2,b.get());
		#else
		operator()(shape,b->bound,b->state->se2,b.get());
		#endif
		if(!b->bound) return; // the functor did not create new bound
		b->bound->refPos=b->state->pos;
		b->bound->lastUpdateIter=scene->iter;
		const Real& sweepLength = b->bound->sweepLength;
		if(sweepLength>0){
			Aabb* aabb=SUDODEM_CAST<Aabb*>(b->bound.get());
			aabb->min-=Vector2r(sweepLength,sweepLength);
			aabb->max+=Vector2r(sweepLength,sweepLength);
		}
	}


/********************************************************************
                      IGeomDispatcher
*********************************************************************/

CREATE_LOGGER(IGeomDispatcher);

shared_ptr<Interaction> IGeomDispatcher::explicitAction(const shared_ptr<Body>& b1, const shared_ptr<Body>& b2, bool force){
	scene=Omega::instance().getScene().get(); // to make sure if called from outside of the loop
	Vector2i cellDist=Vector2i::Zero();
	if(scene->isPeriodic) {
		//throw logic_error("IGeomDispatcher::explicitAction does not support periodic boundary conditions (O.periodic==True)");
		for(int i=0; i<2; i++) cellDist[i]=-(int)((b2->state->pos[i]-b1->state->pos[i])/scene->cell->getSize()[i]+.5);
	}
	Vector2r shift2=scene->cell->hSize*cellDist.cast<Real>();
	updateScenePtr();
	if(force){
		assert(b1->shape && b2->shape);
		shared_ptr<Interaction> I(new Interaction(b1->getId(),b2->getId()));
		I->cellDist=cellDist;
		// FIXME: this code is more or less duplicated from InteractionLoop :-(
		bool swap=false;
		I->functorCache.geom=getFunctor2D(b1->shape,b2->shape,swap);
		if(!I->functorCache.geom) throw invalid_argument("IGeomDispatcher::explicitAction could not dispatch for given types ("+b1->shape->getClassName()+","+b2->shape->getClassName()+").");
		if(swap){I->swapOrder();}
		const shared_ptr<Body>& b1=Body::byId(I->getId1(),scene);
		const shared_ptr<Body>& b2=Body::byId(I->getId2(),scene);
		bool succ=I->functorCache.geom->go(b1->shape,b2->shape,*b1->state,*b2->state,shift2,/*force*/true,I);
		if(!succ) throw logic_error("Functor "+I->functorCache.geom->getClassName()+"::go returned false, even if asked to force IGeom creation. Please report bug.");
		return I;
	} else {
		shared_ptr<Interaction> I(new Interaction(b1->getId(),b2->getId()));
		I->cellDist=cellDist;
		b1->shape && b2->shape && operator()(b1->shape,b2->shape,*b1->state,*b2->state,shift2,/*force*/ false,I);
		return I;
	}
}

void IGeomDispatcher::action(){
	updateScenePtr();

	shared_ptr<BodyContainer>& bodies = scene->bodies;
	const bool isPeriodic(scene->isPeriodic);
	Matrix2r cellHsize; if(isPeriodic) cellHsize=scene->cell->hSize;
	bool removeUnseenIntrs=(scene->interactions->iterColliderLastRun>=0 && scene->interactions->iterColliderLastRun==scene->iter);
	#ifdef SUDODEM_OPENMP
		const long size=scene->interactions->size();
		#pragma omp parallel for
		for(long i=0; i<size; i++){
			const shared_ptr<Interaction>& I=(*scene->interactions)[i];
	#else
		FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
	#endif
			if(removeUnseenIntrs && !I->isReal() && I->iterLastSeen<scene->iter) {
				scene->interactions->requestErase(I);
				continue;
			}
			const shared_ptr<Body>& b1=(*bodies)[I->getId1()];
			const shared_ptr<Body>& b2=(*bodies)[I->getId2()];

			if(!b1 || !b2){ LOG_DEBUG("Body #"<<(b1?I->getId2():I->getId1())<<" vanished, erasing intr #"<<I->getId1()<<"+#"<<I->getId2()<<"!"); scene->interactions->requestErase(I); continue; }

			bool wasReal=I->isReal();
			if (!b1->shape || !b2->shape) { assert(!wasReal); continue; } // some bodies do not have shape
			bool geomCreated;
			if(!isPeriodic){
				geomCreated=operator()(b1->shape, b2->shape, *b1->state, *b2->state, Vector2r::Zero(), /*force*/ false, I);
			} else{
				Vector2r shift2=cellHsize*I->cellDist.cast<Real>();
				geomCreated=operator()(b1->shape, b2->shape, *b1->state, *b2->state, shift2, /*force*/ false, I);
			}
			// reset && erase interaction that existed but now has no geometry anymore
			if(wasReal && !geomCreated){ scene->interactions->requestErase(I); }
	}
}

/********************************************************************
                      IPhysDispatcher
*********************************************************************/


void IPhysDispatcher::explicitAction(shared_ptr<Material>& pp1, shared_ptr<Material>& pp2, shared_ptr<Interaction>& I){
	updateScenePtr();
	if(!I->geom) throw invalid_argument(string(__FILE__)+": explicitAction received interaction without geom.");
	if(!I->functorCache.phys){
		bool dummy;
		I->functorCache.phys=getFunctor2D(pp1,pp2,dummy);
		if(!I->functorCache.phys) throw invalid_argument("IPhysDispatcher::explicitAction did not find a suitable dispatch for types "+pp1->getClassName()+" and "+pp2->getClassName());
		I->functorCache.phys->go(pp1,pp2,I);
	}
}

void IPhysDispatcher::action()
{
	updateScenePtr();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	#ifdef SUDODEM_OPENMP
		const long size=scene->interactions->size();
		#pragma omp parallel for
		for(long i=0; i<size; i++){
			const shared_ptr<Interaction>& interaction=(*scene->interactions)[i];
	#else
		FOREACH(const shared_ptr<Interaction>& interaction, *scene->interactions){
	#endif
			if(interaction->geom){
				shared_ptr<Body>& b1 = (*bodies)[interaction->getId1()];
				shared_ptr<Body>& b2 = (*bodies)[interaction->getId2()];
				bool hadPhys=(interaction->phys.get() != 0);
				operator()(b1->material, b2->material, interaction);
				assert(interaction->phys);
				if(!hadPhys) interaction->iterMadeReal=scene->iter;
			}
		}
}


/********************************************************************
                      LawDispatcher
*********************************************************************/

CREATE_LOGGER(LawDispatcher);
void LawDispatcher::action(){
	updateScenePtr();
	#ifdef SUDODEM_OPENMP
		const long size=scene->interactions->size();
		#pragma omp parallel for
		for(long i=0; i<size; i++){
			const shared_ptr<Interaction>& I=(*scene->interactions)[i];
	#else
		FOREACH(shared_ptr<Interaction> I, *scene->interactions){
	#endif
		if(I->isReal()){
			assert(I->geom); assert(I->phys);
			operator()(I->geom,I->phys,I.get());
			if(!I->isReal() && I->isFresh(scene)) LOG_ERROR("Law functor deleted interaction that was just created. Please report bug: either this message is spurious, or the functor (or something else) is buggy.");
		}
	}
}
