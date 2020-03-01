// 2009 © Václav Šmilauer <eudoxos@arcig.cz>

#include"GeneralIntegratorInsertionSortCollider.hpp"
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/InteractionContainer.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>
#include<sudodem/pkg/dem/Integrator.hpp>
#include<sudodem/pkg/common/Sphere.hpp>

#include<algorithm>
#include<vector>
#include<boost/static_assert.hpp>

SUDODEM_PLUGIN((GeneralIntegratorInsertionSortCollider))
CREATE_LOGGER(GeneralIntegratorInsertionSortCollider);

// STRIDE
bool GeneralIntegratorInsertionSortCollider::isActivated(){
		// activated if number of bodies changes (hence need to refresh collision information)
		// or the time of scheduled run already came, or we were never scheduled yet
		if(!strideActive) return true;
		if(!integrator) return true;
		if(fastestBodyMaxDist<0){fastestBodyMaxDist=0; return true;}
		fastestBodyMaxDist=integrator->maxVelocitySq;
		if(fastestBodyMaxDist>=1 || fastestBodyMaxDist==0) return true;
		if((size_t)BB[0].size!=2*scene->bodies->size()) return true;
		if(scene->interactions->dirty) return true;
		if(scene->doSort) { scene->doSort=false; return true; }
		return false;
	}

void GeneralIntegratorInsertionSortCollider::action(){
	#ifdef ISC_TIMING
		timingDeltas->start();
	#endif

	long nBodies=(long)scene->bodies->size();
	InteractionContainer* interactions=scene->interactions.get();
	scene->interactions->iterColliderLastRun=-1;

	// periodicity changed, force reinit
	if(scene->isPeriodic != periodic){
		for(int i=0; i<3; i++) BB[i].vec.clear();
		periodic=scene->isPeriodic;
	}
	// pre-conditions
		// adjust storage size
		bool doInitSort=false;
		if (doSort) {
			doInitSort=true;
			doSort=false;
		}
		if(BB[0].size!=2*nBodies){
			long BBsize=BB[0].size;
			LOG_DEBUG("Resize bounds containers from "<<BBsize<<" to "<<nBodies*2<<", will std::sort.");
			// bodies deleted; clear the container completely, and do as if all bodies were added (rather slow…)
			// future possibility: insertion sort with such operator that deleted bodies would all go to the end, then just trim bounds
			if(2*nBodies<BBsize){ for(int i=0; i<3; i++) BB[i].vec.clear(); }
			// more than 100 bodies was added, do initial sort again
			// maybe: should rather depend on ratio of added bodies to those already present...?
			if(2*nBodies-BBsize>200 || BBsize==0) doInitSort=true;
			assert((BBsize%2)==0);
			for(int i=0; i<3; i++){
				BB[i].vec.reserve(2*nBodies);
				// add lower and upper bounds; coord is not important, will be updated from bb shortly
				for(long id=BBsize/2; id<nBodies; id++){ BB[i].vec.push_back(Bounds(0,id,/*isMin=*/true)); BB[i].vec.push_back(Bounds(0,id,/*isMin=*/false)); }
				BB[i].size=BB[i].vec.size();
			}
		}
		if(minima.size()!=(size_t)3*nBodies){ minima.resize(3*nBodies); maxima.resize(3*nBodies); }
		assert((size_t)BB[0].size==2*scene->bodies->size());

		// update periodicity
		assert(BB[0].axis==0); assert(BB[1].axis==1); assert(BB[2].axis==2);
		if(periodic) for(int i=0; i<3; i++) BB[i].updatePeriodicity(scene);

		// compatibility block, can be removed later
		findBoundDispatcherInEnginesIfNoFunctorsAndWarn();

		if(verletDist<0){
			Real minR=std::numeric_limits<Real>::infinity();
			FOREACH(const shared_ptr<Body>& b, *scene->bodies){
				if(!b || !b->shape) continue;
				Sphere* s=dynamic_cast<Sphere*>(b->shape.get());
				if(!s) continue;
				minR=min(s->radius,minR);
			}
			if (isinf(minR)) LOG_ERROR("verletDist is set to 0 because no spheres were found. It will result in suboptimal performances, consider setting a positive verletDist in your script.");
			// if no spheres, disable stride
			verletDist=isinf(minR) ? 0 : std::abs(verletDist)*minR;
		}

		// update bounds via boundDispatcher
		boundDispatcher->scene=scene;
		boundDispatcher->sweepDist=verletDist;
		boundDispatcher->minSweepDistFactor=minSweepDistFactor;
		boundDispatcher->targetInterv=targetInterv;
		boundDispatcher->updatingDispFactor=updatingDispFactor;
		boundDispatcher->action();

		// if interactions are dirty, force reinitialization
		if(scene->interactions->dirty){
			doInitSort=true;
			scene->interactions->dirty=false;
		}

		// STRIDE
		if(verletDist>0){
			// get the Integrator, to ask for the maximum velocity value
			if(!integrator){
 				FOREACH(shared_ptr<Engine>& e, scene->engines){ integrator=SUDODEM_PTR_DYN_CAST<Integrator>(e); if(integrator) break; }
				if(!integrator){ throw runtime_error("InsertionSortCollider.verletDist>0, but unable to locate any Integrator within O.engines."); }
			}
		}
	ISC_CHECKPOINT("init");

		// STRIDE
			// get us ready for strides, if they were deactivated
			if(!strideActive && verletDist>0 && integrator->maxVelocitySq>=0){ // maxVelocitySq is a really computed value
				strideActive=true;
			}
			if(strideActive){
				assert(verletDist>0);
				assert(strideActive); assert(integrator->maxVelocitySq>=0);
					integrator->updatingDispFactor=updatingDispFactor;
			} else { /* !strideActive */
				boundDispatcher->sweepDist=0;
			}

	ISC_CHECKPOINT("bound");

	// copy bounds along given axis into our arrays
		for(long i=0; i<2*nBodies; i++){
			for(int j=0; j<3; j++){
				VecBounds& BBj=BB[j];
				const Body::id_t id=BBj[i].id;
				const shared_ptr<Body>& b=Body::byId(id,scene);
				if(b){
					const shared_ptr<Bound>& bv=b->bound;
					// coordinate is min/max if has bounding volume, otherwise both are the position. Add periodic shift so that we are inside the cell
					// watch out for the parentheses around ?: within ?: (there was unwanted conversion of the Reals to bools!)

					BBj[i].coord=((BBj[i].flags.hasBB=((bool)bv)) ? (BBj[i].flags.isMin ? bv->min[j] : bv->max[j]) : (b->state->pos[j])) - (periodic ? BBj.cellDim*BBj[i].period : 0.);

				} else { BBj[i].flags.hasBB=false; /* for vanished body, keep the coordinate as-is, to minimize inversions. */ }
				// if initializing periodic, shift coords & record the period into BBj[i].period
				if(doInitSort && periodic) {
					BBj[i].coord=cellWrap(BBj[i].coord,0,BBj.cellDim,BBj[i].period);
				}
			}
		}
	// for each body, copy its minima and maxima, for quick checks of overlaps later
	BOOST_STATIC_ASSERT(sizeof(Vector3r)==3*sizeof(Real));
	for(Body::id_t id=0; id<nBodies; id++){
		const shared_ptr<Body>& b=Body::byId(id,scene);
		if(b){
			const shared_ptr<Bound>& bv=b->bound;
			if(bv) { memcpy(&minima[3*id],&bv->min,3*sizeof(Real)); memcpy(&maxima[3*id],&bv->max,3*sizeof(Real)); } // ⇐ faster than 6 assignments
			else{ const Vector3r& pos=b->state->pos; memcpy(&minima[3*id],&pos,3*sizeof(Real)); memcpy(&maxima[3*id],&pos,3*sizeof(Real)); }
		} else { memset(&minima[3*id],0,3*sizeof(Real)); memset(&maxima[3*id],0,3*sizeof(Real)); }
	}

	ISC_CHECKPOINT("copy");

	// process interactions that the constitutive law asked to be erased
// 	interactions->erasePending(*this,scene);
	interactions->conditionalyEraseNonReal(*this,scene);

	ISC_CHECKPOINT("erase");

	// sort
		// the regular case
		if(!doInitSort && !sortThenCollide){
			/* each inversion in insertionSort calls handleBoundInversion, which in turns may add/remove interaction */
			if(!periodic) for(int i=0; i<3; i++) insertionSort(BB[i],interactions,scene);
			else for(int i=0; i<3; i++) insertionSortPeri(BB[i],interactions,scene);
		}
		// create initial interactions (much slower)
		else {
			if(doInitSort){
				// the initial sort is in independent in 3 dimensions, may be run in parallel; it seems that there is no time gain running in parallel, though
				// important to reset loInx for periodic simulation (!!)
				for(int i=0; i<3; i++) { BB[i].loIdx=0; std::sort(BB[i].vec.begin(),BB[i].vec.end()); }
				numReinit++;
			} else { // sortThenCollide
				if(!periodic) for(int i=0; i<3; i++) insertionSort(BB[i],interactions,scene,false);
				else for(int i=0; i<3; i++) insertionSortPeri(BB[i],interactions,scene,false);
			}
			// traverse the container along requested axis
			assert(sortAxis==0 || sortAxis==1 || sortAxis==2);
			VecBounds& V=BB[sortAxis];
			// go through potential aabb collisions, create interactions as necessary
			if(!periodic){
				for(long i=0; i<2*nBodies; i++){
					// start from the lower bound (i.e. skipping upper bounds)
					// skip bodies without bbox, because they don't collide
					if(!(V[i].flags.isMin && V[i].flags.hasBB)) continue;
					const Body::id_t& iid=V[i].id;
					// go up until we meet the upper bound
					for(long j=i+1; /* handle case 2. of swapped min/max */ j<2*nBodies && V[j].id!=iid; j++){
						const Body::id_t& jid=V[j].id;
						// take 2 of the same condition (only handle collision [min_i..max_i]+min_j, not [min_i..max_i]+min_i (symmetric)
						if(!V[j].flags.isMin) continue;
						/* abuse the same function here; since it does spatial overlap check first, it is OK to use it */
						handleBoundInversion(iid,jid,interactions,scene);
						assert(j<2*nBodies-1);
					}
				}
			} else { // periodic case: see comments above
				for(long i=0; i<2*nBodies; i++){
					if(!(V[i].flags.isMin && V[i].flags.hasBB)) continue;
					const Body::id_t& iid=V[i].id;
					long cnt=0;
					// we might wrap over the periodic boundary here; that's why the condition is different from the aperiodic case
					for(long j=V.norm(i+1); V[j].id!=iid; j=V.norm(j+1)){
						const Body::id_t& jid=V[j].id;
						if(!V[j].flags.isMin) continue;
						handleBoundInversionPeri(iid,jid,interactions,scene);
						if(cnt++>2*(long)nBodies){ LOG_FATAL("Uninterrupted loop in the initial sort?"); throw std::logic_error("loop??"); }
					}
				}
			}
		}
	ISC_CHECKPOINT("sort&collide");
}
