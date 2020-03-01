#include"InteractionLoop.hpp"

SUDODEM_PLUGIN((InteractionLoop));
CREATE_LOGGER(InteractionLoop);

void InteractionLoop::pyHandleCustomCtorArgs(boost::python::tuple& t, boost::python::dict& d){
	if(boost::python::len(t)==0) return; // nothing to do
	if(boost::python::len(t)!=3) throw invalid_argument("Exactly 3 lists of functors must be given");
	// parse custom arguments (3 lists) and do in-place modification of args
	typedef std::vector<shared_ptr<IGeomFunctor> > vecGeom;
	typedef std::vector<shared_ptr<IPhysFunctor> > vecPhys;
	typedef std::vector<shared_ptr<LawFunctor> > vecLaw;
	vecGeom vg=boost::python::extract<vecGeom>(t[0])();
	vecPhys vp=boost::python::extract<vecPhys>(t[1])();
	vecLaw vl=boost::python::extract<vecLaw>(t[2])();
	FOREACH(shared_ptr<IGeomFunctor> gf, vg) this->geomDispatcher->add(gf);
	FOREACH(shared_ptr<IPhysFunctor> pf, vp) this->physDispatcher->add(pf);
	FOREACH(shared_ptr<LawFunctor> cf, vl) this->lawDispatcher->add(cf);
	t=boost::python::tuple(); // empty the args; not sure if this is OK, as there is some refcounting in raw_constructor code
}


void InteractionLoop::action(){
	// update Scene* of the dispatchers
	geomDispatcher->scene=physDispatcher->scene=lawDispatcher->scene=scene;
	// ask dispatchers to update Scene* of their functors
	geomDispatcher->updateScenePtr(); physDispatcher->updateScenePtr(); lawDispatcher->updateScenePtr();

	// call Ig2Functor::preStep
	FOREACH(const shared_ptr<IGeomFunctor>& ig2, geomDispatcher->functors) ig2->preStep();
	// call LawFunctor::preStep
	FOREACH(const shared_ptr<LawFunctor>& law2, lawDispatcher->functors) law2->preStep();

	/*
		initialize callbacks; they return pointer (used only in this timestep) to the function to be called
		returning NULL deactivates the callback in this timestep
	*/
	// pair of callback object and pointer to the function to be called
	vector<IntrCallback::FuncPtr> callbackPtrs;
	FOREACH(const shared_ptr<IntrCallback> cb, callbacks){
		cb->scene=scene;
		callbackPtrs.push_back(cb->stepInit());
	}
	assert(callbackPtrs.size()==callbacks.size());
	size_t callbacksSize=callbacks.size();

	// cache transformed cell size
	Matrix2r cellHsize; if(scene->isPeriodic) cellHsize=scene->cell->hSize;

	// force removal of interactions that were not encountered by the collider
	// (only for some kinds of colliders; see comment for InteractionContainer::iterColliderLastRun)
	bool removeUnseenIntrs=(scene->interactions->iterColliderLastRun>=0 && scene->interactions->iterColliderLastRun==scene->iter);

	#ifdef SUDODEM_OPENMP
	const long size=scene->interactions->size();
	#pragma omp parallel for schedule(guided) num_threads(ompThreads>0 ? min(ompThreads,omp_get_max_threads()) : omp_get_max_threads())
	for(long i=0; i<size; i++){
		const shared_ptr<Interaction>& I=(*scene->interactions)[i];
	#else
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
	#endif
		// keep the following newline, my (edx) preprocessor outputs garbage code otherwise!

		if(removeUnseenIntrs && !I->isReal() && I->iterLastSeen<scene->iter) {
			eraseAfterLoop(I->getId1(),I->getId2());
			continue;
		}

		const shared_ptr<Body>& b1_=Body::byId(I->getId1(),scene);
		const shared_ptr<Body>& b2_=Body::byId(I->getId2(),scene);

		if(!b1_ || !b2_){ LOG_DEBUG("Body #"<<(b1_?I->getId2():I->getId1())<<" vanished, erasing intr #"<<I->getId1()<<"+#"<<I->getId2()<<"!"); scene->interactions->requestErase(I); continue; }

    // Skip interaction with clumps
    if (b1_->isClump() || b2_->isClump()) { continue; }
		// we know there is no geometry functor already, take the short path
		if(!I->functorCache.geomExists) { assert(!I->isReal()); continue; }
		// no interaction geometry for either of bodies; no interaction possible
		if(!b1_->shape || !b2_->shape) { assert(!I->isReal()); continue; }

		bool swap=false;
		// IGeomDispatcher
		if(!I->functorCache.geom){
			I->functorCache.geom=geomDispatcher->getFunctor2D(b1_->shape,b2_->shape,swap);
			// returns NULL ptr if no functor exists; remember that and shortcut
			if(!I->functorCache.geom) {I->functorCache.geomExists=false; continue; }

		}
		// arguments for the geom functor are in the reverse order (dispatcher would normally call goReverse).
		// we don't remember the fact that is reverse, so we swap bodies within the interaction
		// and can call go in all cases
		if(swap){I->swapOrder();}
		// body pointers must be updated, in case we swapped
		const shared_ptr<Body>& b1=swap?b2_:b1_;
		const shared_ptr<Body>& b2=swap?b1_:b2_;

		assert(I->functorCache.geom);
		bool wasReal=I->isReal();
		bool geomCreated;
		if(!scene->isPeriodic){
			geomCreated=I->functorCache.geom->go(b1->shape,b2->shape, *b1->state, *b2->state, Vector2r::Zero(), /*force*/false, I);
		} else { // handle periodicity
			Vector2r shift2=cellHsize*I->cellDist.cast<Real>();
			// in sheared cell, apply shear on the mutual position as well
			geomCreated=I->functorCache.geom->go(b1->shape,b2->shape,*b1->state,*b2->state,shift2,/*force*/false,I);
		}

		if(!geomCreated){
			if(wasReal) LOG_WARN("IGeomFunctor returned false on existing interaction!");
			if(wasReal) scene->interactions->requestErase(I); // fully created interaction without geometry is reset and perhaps erased in the next step
			continue; // in any case don't care about this one anymore
		}

		// IPhysDispatcher
		if(!I->functorCache.phys){
			I->functorCache.phys=physDispatcher->getFunctor2D(b1->material,b2->material,swap);
			assert(!swap); // InteractionPhysicsEngineUnits are symmetric
		}
		//assert(I->functorCache.phys);
		if(!I->functorCache.phys){
			throw std::runtime_error("Undefined or ambiguous IPhys dispatch for types "+b1->material->getClassName()+" and "+b2->material->getClassName()+".");
		}
		I->functorCache.phys->go(b1->material,b2->material,I);
		
		assert(I->phys);
		if(!wasReal) I->iterMadeReal=scene->iter; // mark the interaction as created right now

		// LawDispatcher
		// populating constLaw cache must be done after geom and physics dispatchers have been called, since otherwise the interaction
		// would not have geom and phys yet.
		if(!I->functorCache.constLaw){
			I->functorCache.constLaw=lawDispatcher->getFunctor2D(I->geom,I->phys,swap);
			if(!I->functorCache.constLaw){
				LOG_FATAL("None of given Law2 functors can handle interaction #"<<I->getId1()<<"+"<<I->getId2()<<", types geom:"<<I->geom->getClassName()<<"="<<I->geom->getClassIndex()<<" and phys:"<<I->phys->getClassName()<<"="<<I->phys->getClassIndex()<<" (LawDispatcher::getFunctor2D returned empty functor)");
				//abort();
				exit(1);
			}
			assert(!swap); // reverse call would make no sense, as the arguments are of different types
		}
		assert(I->functorCache.constLaw);
		//If the functor return false, the interaction is reset
		if (!I->functorCache.constLaw->go(I->geom,I->phys,I.get())) scene->interactions->requestErase(I);

		// process callbacks for this interaction
		if(!I->isReal()) continue; // it is possible that Law2_ functor called requestErase, hence this check
		for(size_t i=0; i<callbacksSize; i++){
			if(callbackPtrs[i]!=NULL) (*(callbackPtrs[i]))(callbacks[i].get(),I.get());
		}
	}
}
