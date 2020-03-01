// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
// 2013 © Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>

#include"InsertionSortCollider.hpp"
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/InteractionContainer.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>
#include<sudodem/pkg/dem/NewtonIntegrator.hpp>
#include<sudodem/pkg/common/Disk.hpp>

#include<boost/static_assert.hpp>
#ifdef SUDODEM_OPENMP
  #include<omp.h>
#endif

//#include<sudodem/lib/serialization/ObjectIO.hpp>//output state for debuging

SUDODEM_PLUGIN((InsertionSortCollider))
CREATE_LOGGER(InsertionSortCollider);


// called by the insertion sort if 2 bodies swapped their bounds in such a way that a new overlap may appear
void InsertionSortCollider::handleBoundInversion(Body::id_t id1, Body::id_t id2, InteractionContainer* interactions, Scene*){
	assert(!periodic);
	assert(id1!=id2);
	if (spatialOverlap(id1,id2) && Collider::mayCollide(Body::byId(id1,scene).get(),Body::byId(id2,scene).get()) && !interactions->found(id1,id2))
		interactions->insert(shared_ptr<Interaction>(new Interaction(id1,id2)));
}

void InsertionSortCollider::insertionSort(VecBounds& v, InteractionContainer* interactions, Scene*, bool doCollide){
	assert(!periodic);
	assert(v.size==(long)v.vec.size());
	for(long i=1; i<v.size; i++){
		const Bounds viInit=v[i]; long j=i-1; /* cache hasBB; otherwise 1% overall performance hit */ const bool viInitBB=viInit.flags.hasBB;
		const bool isMin=viInit.flags.isMin;

		while(j>=0 && v[j]>viInit){
			v[j+1]=v[j];
			// no collisions without bounding boxes
			// also, do not collide body with itself; it sometimes happens for facets aligned perpendicular to an axis, for reasons that are not very clear
			// see https://bugs.launchpad.net/sudodem/+bug/669095
			// skip bounds with same isMin flags, since inversion doesn't imply anything in that case
			if(isMin && !v[j].flags.isMin && doCollide && viInitBB && v[j].flags.hasBB && (viInit.id!=v[j].id)) {
				/*if (isMin)*/ handleBoundInversion(viInit.id,v[j].id,interactions,scene);
// 				else handleBoundSplit(viInit.id,v[j].id,interactions,scene);
			}
			j--;
		}
		v[j+1]=viInit;
	}
}


//Periodic version, only for non-periodic case at the moment (feel free to implement for the periodic case...)
void InsertionSortCollider::insertionSortParallel(VecBounds& v, InteractionContainer* interactions, Scene*, bool doCollide){
#ifdef SUDODEM_OPENMP
	assert(!periodic);
	assert(v.size==(long)v.vec.size());
	if (ompThreads<=1) return insertionSort(v,interactions, scene, doCollide);

	Real chunksVerlet = 4*verletDist;//is 2* the theoretical requirement?
	if (chunksVerlet<=0) {LOG_ERROR("Parallel insertion sort needs verletDist>0");}

	///chunks defines subsets of the bounds lists, we make sure they are not too small wrt. verlet dist.
	std::vector<unsigned> chunks;
	unsigned nChunks = ompThreads;
	unsigned chunkSize = unsigned(v.size/nChunks)+1;
	for(unsigned n=0; n<nChunks;n++) chunks.push_back(n*chunkSize); chunks.push_back(v.size);
	while (nChunks>1){
		bool changeChunks=false;
		for(unsigned n=1; n<nChunks;n++) if (chunksVerlet>(v[chunks[n]].coord-v[chunks[n-1]].coord)) changeChunks=true;
		if (!changeChunks) break;
		nChunks--; chunkSize = unsigned(v.size/nChunks)+1; chunks.clear();
		for(unsigned n=0; n<nChunks;n++) chunks.push_back(n*chunkSize); chunks.push_back(v.size);
	}
	static unsigned warnOnce=0;
	if (nChunks<unsigned(ompThreads) && !warnOnce++) LOG_WARN("Parallel insertion: only "<<nChunks <<" thread(s) used. The number of bodies is probably too small for allowing more threads. The contact detection should succeed but not all available threads are used.");

	///Define per-thread containers bufferizing the actual insertion of new interactions, since inserting is not thread-safe
	std::vector<std::vector<std::pair<Body::id_t,Body::id_t> > > newInteractions;
	newInteractions.resize(ompThreads,std::vector<std::pair<Body::id_t,Body::id_t> >());
	for (int kk=0;  kk<ompThreads; kk++) newInteractions[kk].reserve(100);

	/// First sort, independant in each chunk
	#pragma omp parallel for schedule(dynamic,1) num_threads(ompThreads>0 ? min(ompThreads,omp_get_max_threads()) : omp_get_max_threads())
	for (unsigned k=0; k<nChunks;k++) {
		int threadNum = omp_get_thread_num();
		for(long i=chunks[k]+1; i<chunks[k+1]; i++){
			const Bounds viInit=v[i]; long j=i-1; const bool viInitBB=viInit.flags.hasBB;
			const bool isMin=viInit.flags.isMin;
			while(j>=chunks[k] && v[j]>viInit){
				v[j+1]=v[j];
				if(isMin && !v[j].flags.isMin && doCollide && viInitBB && v[j].flags.hasBB && (viInit.id!=v[j].id)) {
					const Body::id_t& id1 = v[j].id; const Body::id_t& id2 = viInit.id;
					if (spatialOverlap(id1,id2) && Collider::mayCollide(Body::byId(id1,scene).get(),Body::byId(id2,scene).get()) && !interactions->found(id1,id2))
						newInteractions[threadNum].push_back(std::pair<Body::id_t,Body::id_t>(v[j].id,viInit.id));
				}
				j--;
			}
			v[j+1]=viInit;
		}
	}

	///In the second sort, the chunks are connected consistently.
	///If sorting requires to move a bound past half-chunk, the algorithm is not thread safe,
	/// if it happens we roll-back and run the 1-thread sort + send warning
	bool parallelFailed=false;
	#pragma omp parallel for schedule(dynamic,1) num_threads(ompThreads>0 ? min(ompThreads,omp_get_max_threads()) : omp_get_max_threads())
	for (unsigned k=1; k<nChunks;k++) {

		int threadNum = omp_get_thread_num();
		long i=chunks[k];
		for(; v[i]<v[i-1]; i++){
			const Bounds viInit=v[i]; long j=i-1; /* cache hasBB; otherwise 1% overall performance hit */ const bool viInitBB=viInit.flags.hasBB;
			const bool isMin=viInit.flags.isMin;

			while(j>=chunks[k-1] && viInit<v[j]){
				v[j+1]=v[j];
				if(isMin && !v[j].flags.isMin && doCollide && viInitBB && v[j].flags.hasBB && (viInit.id!=v[j].id)) {
					const Body::id_t& id1 = v[j].id; const Body::id_t& id2 = viInit.id;
					//FIXME: do we need the check with found(id1,id2) here? It is checked again below...
					if (spatialOverlap(id1,id2) && Collider::mayCollide(Body::byId(id1,scene).get(),Body::byId(id2,scene).get()) && !interactions->found(id1,id2))
						newInteractions[threadNum].push_back(std::pair<Body::id_t,Body::id_t>(v[j].id,viInit.id));}
				j--;
			}
			v[j+1]=viInit;
			if (j<=long(chunks[k]-chunkSize*0.5)) {
				LOG_WARN("parallel sort not guaranteed to succeed; in chunk "<<k<<" of "<<nChunks<< ", bound descending past half-chunk. Parallel colliding aborted, starting again in single thread. Consider turning collider.ompThreads=1 for not wasting CPU time.");
				parallelFailed=true;}
		}
		if (i>=long(chunks[k]+chunkSize*0.5)) {
			LOG_ERROR("parallel sort not guaranteed to succeed; in chunk "<<k+1<<" of "<<nChunks<< ", bound advancing past half-chunk. Consider turning collider.ompThreads=1 for not wasting CPU time.");
			parallelFailed=true;}
	}
	/// Check again, just to be sure...
	for (unsigned k=1; k<nChunks;k++) if (v[chunks[k]]<v[chunks[k]-1]) {
		LOG_ERROR("Parallel colliding failed, starting again in single thread. Consider turning collider.ompThreads=1 for not wasting CPU time.");
		parallelFailed=true;}

	if (parallelFailed) return insertionSort(v,interactions, scene, doCollide);

	/// Now insert interactions sequentially
	for (int n=0;n<ompThreads;n++) for (size_t k=0, kend=newInteractions[n].size();k<kend;k++) if (!interactions->found(newInteractions[n][k].first,newInteractions[n][k].second)) interactions->insert(shared_ptr<Interaction>(new Interaction(newInteractions[n][k].first,newInteractions[n][k].second)));
#endif
}


vector<Body::id_t> InsertionSortCollider::probeBoundingVolume(const Bound& bv){
	if(periodic){ throw invalid_argument("InsertionSortCollider::probeBoundingVolume: handling periodic boundary not implemented."); }
	vector<Body::id_t> ret;
	for( vector<Bounds>::iterator
			it=BB[0].vec.begin(),et=BB[0].vec.end(); it < et; ++it)
	{
		if (it->coord > bv.max[0]) break;
		if (!it->flags.isMin || !it->flags.hasBB) continue;
		int offset = 2*it->id;
		const shared_ptr<Body>& b=Body::byId(it->id,scene);
		if(!b || !b->bound) continue;
		const Real& sweepLength = b->bound->sweepLength;
		Vector2r disp = b->state->pos - b->bound->refPos;
		if (!(maxima[offset]-sweepLength+disp[0] < bv.min[0] ||
			minima[offset]+sweepLength+disp[0] > bv.max[0] ||
			minima[offset+1]+sweepLength+disp[1] > bv.max[1] ||
			maxima[offset+1]-sweepLength+disp[1] < bv.min[1]
			))
		{
			ret.push_back(it->id);
		}
	}
	return ret;
}
// STRIDE
	bool InsertionSortCollider::isActivated(){
		// activated if number of bodies changes (hence need to refresh collision information)
		// or the time of scheduled run already came, or we were never scheduled yet
		if(!strideActive) return true;
		if(!newton) return true;
		if(fastestBodyMaxDist<0){fastestBodyMaxDist=0; return true;}
		fastestBodyMaxDist=newton->maxVelocitySq;
		if(fastestBodyMaxDist>=1 || fastestBodyMaxDist==0) return true;
		if((size_t)BB[0].size!=2*scene->bodies->size()) return true;
		if(scene->interactions->dirty) return true;
		if(scene->doSort) { scene->doSort=false; return true; }
		return false;
	}

void InsertionSortCollider::action(){
	#ifdef ISC_TIMING
		timingDeltas->start();
	#endif

	long nBodies=(long)scene->bodies->size();
	InteractionContainer* interactions=scene->interactions.get();
	scene->interactions->iterColliderLastRun=-1;
	#ifdef SUDODEM_OPENMP
	if (ompThreads<=0) ompThreads = omp_get_max_threads();
	#endif
	// periodicity changed, force reinit
	if(scene->isPeriodic != periodic){
		for(int i=0; i<2; i++) BB[i].vec.clear();
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
			if(2*nBodies<BBsize){ for(int i=0; i<2; i++) BB[i].vec.clear(); }
			// more than 100 bodies was added, do initial sort again
			// maybe: should rather depend on ratio of added bodies to those already present...?
			if(2*nBodies-BBsize>200 || BBsize==0) doInitSort=true;
			assert((BBsize%2)==0);
			for(int i=0; i<2; i++){
				BB[i].vec.reserve(2*nBodies);
				// add lower and upper bounds; coord is not important, will be updated from bb shortly
				for(long id=BBsize/2; id<nBodies; id++){ BB[i].vec.push_back(Bounds(0,id,/*isMin=*/true)); BB[i].vec.push_back(Bounds(0,id,/*isMin=*/false)); }
				BB[i].size=BB[i].vec.size();
			}
		}
		if(minima.size()!=(size_t)2*nBodies){ minima.resize(2*nBodies); maxima.resize(2*nBodies); }
		assert((size_t)BB[0].size==2*scene->bodies->size());
    //Increase the size of force container.
    scene->forces.addMaxId(2*scene->bodies->size());//need?
		// update periodicity
		assert(BB[0].axis==0); assert(BB[1].axis==1);
		if(periodic){
       for(int i=0; i<2; i++) BB[i].updatePeriodicity(scene);
       invSizes=Vector2r(1./scene->cell->getSize()[0],1./scene->cell->getSize()[1]);
    }
		// compatibility block, can be removed later
		//findBoundDispatcherInEnginesIfNoFunctorsAndWarn();

		if(verletDist<0){
			Real minR=std::numeric_limits<Real>::infinity();
			FOREACH(const shared_ptr<Body>& b, *scene->bodies){
				if(!b || !b->shape) continue;
				Disk* s=dynamic_cast<Disk*>(b->shape.get());//fixme for other shapes
				if(!s) continue;
				minR=min(s->radius,minR);
			}
			if (isinf(minR)) LOG_ERROR("verletDist is set to 0 because no disks were found. It will result in suboptimal performances, consider setting a positive verletDist in your script.");
			// if no disks, disable stride
			verletDist=isinf(minR) ? 0 : std::abs(verletDist)*minR;
		}
		// if interactions are dirty, force reinitialization
		if(scene->interactions->dirty){
			doInitSort=true;
			scene->interactions->dirty=false;
		}
    //cout<<"hello:updating collider!"<<endl;
		// update bounds via boundDispatcher
		boundDispatcher->scene=scene;
		boundDispatcher->sweepDist=verletDist;
		boundDispatcher->minSweepDistFactor=minSweepDistFactor;
		boundDispatcher->targetInterv=targetInterv;
		boundDispatcher->updatingDispFactor=updatingDispFactor;
		boundDispatcher->action();
		ISC_CHECKPOINT("boundDispatcher");


		// STRIDE
		if(verletDist>0){
			// get NewtonIntegrator, to ask for the maximum velocity value
			if(!newton){
				FOREACH(shared_ptr<Engine>& e, scene->engines){ newton=SUDODEM_PTR_DYN_CAST<NewtonIntegrator>(e); if(newton) break; }
				if(!newton){ throw runtime_error("InsertionSortCollider.verletDist>0, but unable to locate NewtonIntegrator within O.engines."); }
			}
		}
		// STRIDE
			// get us ready for strides, if they were deactivated
			if(!strideActive && verletDist>0 && newton->maxVelocitySq>=0) strideActive=true;
			if(strideActive){
				assert(verletDist>0);
				assert(strideActive); assert(newton->maxVelocitySq>=0);
				newton->updatingDispFactor=updatingDispFactor;
			} else boundDispatcher->sweepDist=0;

	ISC_CHECKPOINT("bound");

	// copy bounds along given axis into our arrays
	#pragma omp parallel for schedule(guided) num_threads(ompThreads>0 ? min(ompThreads,omp_get_max_threads()) : omp_get_max_threads())
	for(long i=0; i<2*nBodies; i++){
// 		const long cacheIter = scene->iter;
		for(int j=0; j<2; j++){
				VecBounds& BBj=BB[j];
				Bounds& BBji = BBj[i];
				const Body::id_t id=BBji.id;
				const shared_ptr<Body>& b=Body::byId(id,scene);
				if(b){
					const shared_ptr<Bound>& bv=b->bound;
					// coordinate is min/max if has bounding volume, otherwise both are the position. Add periodic shift so that we are inside the cell
					// watch out for the parentheses around ?: within ?: (there was unwanted conversion of the Reals to bools!)
					BBji.coord=((BBji.flags.hasBB=((bool)bv)) ? (BBji.flags.isMin ? bv->min[j] : bv->max[j]) : (b->state->pos[j])) - (periodic ? BBj.cellDim*BBji.period : 0.);
					// if initializing periodic, shift coords & record the period into BBj[i].period
					if(doInitSort && periodic) BBji.coord=cellWrap(BBji.coord,0,BBj.cellDim,BBji.period);
					// for each body, copy its minima and maxima, for quick checks of overlaps later
					//bounds have been all updated when j==0, we can safely copy them here when j==1
					if (BBji.flags.isMin && j==1 &&bv) {
						 memcpy(&minima[2*id],&bv->min,2*sizeof(Real)); memcpy(&maxima[2*id],&bv->max,2*sizeof(Real));
					}
				} else { BBj[i].flags.hasBB=false; /* for vanished body, keep the coordinate as-is, to minimize inversions. */ }
			}
		}

	ISC_CHECKPOINT("copy");

	// remove interactions which have disconnected bounds and are not real (will run parallel if SUDODEM_OPENMP)
	interactions->conditionalyEraseNonReal(*this,scene);

	ISC_CHECKPOINT("erase");

	// sort
		// the regular case
		if(!doInitSort && !sortThenCollide){
			/* each inversion in insertionSort calls may add interaction */
			//1000 bodies is heuristic minimum above which parallel sort is called
			if(!periodic) for(int i=0; i<2; i++) {
			#ifdef SUDODEM_OPENMP
				if (ompThreads<=1 || nBodies<1000) insertionSort(BB[i],interactions,scene);
				else insertionSortParallel(BB[i],interactions,scene);}
			#else
				insertionSort(BB[i],interactions,scene);}
			#endif
			else for(int i=0; i<2; i++) insertionSortPeri(BB[i],interactions,scene);
		}
		// create initial interactions (much slower)
		else {
			if(doInitSort){
				// the initial sort is in independent in 3 dimensions, may be run in parallel; it seems that there is no time gain running in parallel, though
				// important to reset loInx for periodic simulation (!!)
// 				#pragma omp parallel for schedule(dynamic,1) num_threads(min(ompThreads,3))
				for(int i=0; i<2; i++) { BB[i].loIdx=0; std::sort(BB[i].vec.begin(),BB[i].vec.end()); }
				numReinit++;
			} else { // sortThenCollide
				if(!periodic) for(int i=0; i<2; i++) insertionSort(BB[i],interactions,scene,false);
				else for(int i=0; i<2; i++) insertionSortPeri(BB[i],interactions,scene,false);
			}
			// traverse the container along requested axis
			assert(sortAxis==0 || sortAxis==1);
			VecBounds& V=BB[sortAxis];
			// go through potential aabb collisions, create interactions as necessary
			if(!periodic){
			#ifdef SUDODEM_OPENMP
				std::vector<std::vector<std::pair<Body::id_t,Body::id_t> > > newInts;
				newInts.resize(ompThreads,std::vector<std::pair<Body::id_t,Body::id_t> >());
				for (int kk=0;  kk<ompThreads; kk++) newInts[kk].reserve(unsigned(10*nBodies/ompThreads));
				#pragma omp parallel for schedule(guided,200) num_threads(ompThreads)
			#endif
				for(long i=0; i<2*nBodies; i++){
					// start from the lower bound (i.e. skipping upper bounds)
					// skip bodies without bbox, because they don't collide
					if(!(V[i].flags.isMin && V[i].flags.hasBB)) continue;
					const Body::id_t& iid=V[i].id;
					// go up until we meet the upper bound
					for(long j=i+1; /* handle case 2. of swapped min/max */ j<2*nBodies && V[j].id!=iid; j++){
						const Body::id_t& jid=V[j].id;
						// take 2 of the same condition (only handle collision [min_i..max_i]+min_j, not [min_i..max_i]+min_i (symmetric)
						if(!(V[j].flags.isMin && V[j].flags.hasBB)) continue;
						if (spatialOverlap(iid,jid) && Collider::mayCollide(Body::byId(iid,scene).get(),Body::byId(jid,scene).get()) ){
						#ifdef SUDODEM_OPENMP
							unsigned int threadNum = omp_get_thread_num();
							newInts[threadNum].push_back(std::pair<Body::id_t,Body::id_t>(iid,jid));
						#else
							if (!interactions->found(iid,jid))
							interactions->insert(shared_ptr<Interaction>(new Interaction(iid,jid)));
						#endif
						}
					}
				}
				//go through newly created candidates sequentially, duplicates coming from different threads may exist so we check existence with found()
				#ifdef SUDODEM_OPENMP
				for (int n=0;n<ompThreads;n++) for (size_t k=0, kend=newInts[n].size();k<kend;k++)
					if (!interactions->found(newInts[n][k].first,newInts[n][k].second))
						interactions->insert(shared_ptr<Interaction>(new Interaction(newInts[n][k].first,newInts[n][k].second)));
				#endif
			} else { // periodic case: see comments above
				for(long i=0; i<2*nBodies; i++){
					if(!(V[i].flags.isMin && V[i].flags.hasBB)) continue;
					const Body::id_t& iid=V[i].id;
					long cnt=0;
					// we might wrap over the periodic boundary here; that's why the condition is different from the aperiodic case
					for(long j=V.norm(i+1); V[j].id!=iid; j=V.norm(j+1)){
						const Body::id_t& jid=V[j].id;
						if(!(V[j].flags.isMin && V[j].flags.hasBB)) continue;
						handleBoundInversionPeri(iid,jid,interactions,scene);
						if(cnt++>2*(long)nBodies){ LOG_FATAL("Uninterrupted loop in the initial sort?"); throw std::logic_error("loop??"); }
					}
				}
			}
		}
	ISC_CHECKPOINT("sort&collide");
}


// return floating value wrapped between x0 and x1 and saving period number to period
Real InsertionSortCollider::cellWrap(const Real x, const Real x0, const Real x1, int& period){
	Real xNorm=(x-x0)/(x1-x0);
	period=(int)floor(xNorm); // some people say this is very slow; probably optimized by gcc, however (google suggests)
	return x0+(xNorm-period)*(x1-x0);
}

// return coordinate wrapped to 0…x1, relative to x0; don't care about period
Real InsertionSortCollider::cellWrapRel(const Real x, const Real x0, const Real x1){
	Real xNorm=(x-x0)/(x1-x0);
	return (xNorm-floor(xNorm))*(x1-x0);
}

void InsertionSortCollider::insertionSortPeri(VecBounds& v, InteractionContainer* interactions, Scene*, bool doCollide){
	assert(periodic);
	long &loIdx=v.loIdx; const long &size=v.size;
	for(long _i=0; _i<size; _i++){
		const long i=v.norm(_i);
		const long i_1=v.norm(i-1);
		//switch period of (i) if the coord is below the lower edge cooridnate-wise and just above the split
		if(i==loIdx && v[i].coord<0){ v[i].period-=1; v[i].coord+=v.cellDim; loIdx=v.norm(loIdx+1); }
		// coordinate of v[i] used to check inversions
		// if crossing the split, adjust by cellDim;
		// if we get below the loIdx however, the v[i].coord will have been adjusted already, no need to do that here
		const Real iCmpCoord=v[i].coord+(i==loIdx ? v.cellDim : 0);
		// no inversion
		if(v[i_1].coord<=iCmpCoord) continue;
		// vi is the copy that will travel down the list, while other elts go up
		// if will be placed in the list only at the end, to avoid extra copying
		int j=i_1; Bounds vi=v[i];  const bool viHasBB=vi.flags.hasBB;
		const bool isMin=v[i].flags.isMin;
		while(v[j].coord>vi.coord + /* wrap for elt just below split */ (v.norm(j+1)==loIdx ? v.cellDim : 0)){
			long j1=v.norm(j+1);
			// OK, now if many bodies move at the same pace through the cell and at one point, there is inversion,
			// this can happen without any side-effects
			if (false && v[j].coord>2*v.cellDim){
				// this condition is not strictly necessary, but the loop of insertionSort would have to run more times.
				// Since size of particle is required to be < .5*cellDim, this would mean simulation explosion anyway
				LOG_FATAL("Body #"<<v[j].id<<" going faster than 1 cell in one step? Not handled.");
				throw runtime_error(__FILE__ ": body mmoving too fast (skipped 1 cell).");
			}
			Bounds& vNew(v[j1]); // elt at j+1 being overwritten by the one at j and adjusted
			vNew=v[j];
			// inversions close the the split need special care
			if(j==loIdx && vi.coord<0) { vi.period-=1; vi.coord+=v.cellDim; loIdx=v.norm(loIdx+1); }
			else if(j1==loIdx) { vNew.period+=1; vNew.coord-=v.cellDim; loIdx=v.norm(loIdx-1); }
			if(isMin && !v[j].flags.isMin && (doCollide && viHasBB && v[j].flags.hasBB)){
				// see https://bugs.launchpad.net/sudodem/+bug/669095 and similar problem in aperiodic insertionSort
				#if 0
				if(vi.id==vNew.id){
					LOG_FATAL("Inversion of body's #"<<vi.id<<" boundary with its other boundary, "<<v[j].coord<<" meets "<<vi.coord);
					throw runtime_error(__FILE__ ": Body's boundary metting its opposite boundary.");
				}
				#endif
				if((vi.id!=vNew.id)) handleBoundInversionPeri(vi.id,vNew.id,interactions,scene);
			}
			j=v.norm(j-1);
		}
		v[v.norm(j+1)]=vi;
	}
}

// called by the insertion sort if 2 bodies swapped their bounds
void InsertionSortCollider::handleBoundInversionPeri(Body::id_t id1, Body::id_t id2, InteractionContainer* interactions, Scene*){
	assert(periodic);
	Vector2i periods;
	bool overlap=spatialOverlapPeri(id1,id2,scene,periods);
	if (overlap && Collider::mayCollide(Body::byId(id1,scene).get(),Body::byId(id2,scene).get()) && !interactions->found(id1,id2)){
		shared_ptr<Interaction> newI=shared_ptr<Interaction>(new Interaction(id1,id2));
		newI->cellDist=periods;
		interactions->insert(newI);
	}
}

/* Performance hint
	================

	Since this function is called (most the time) from insertionSort,
	we actually already know what is the overlap status in that one dimension, including
	periodicity information; both values could be passed down as a parameters, avoiding 1 of 3 loops here.
	We do some floats math here, so the speedup could noticeable; doesn't concertn the non-periodic variant,
	where it is only plain comparisons taking place.

	In the same way, handleBoundInversion is passed only id1 and id2, but again insertionSort already knows in which sense
	the inversion happens; if the boundaries get apart (min getting up over max), it wouldn't have to check for overlap
	at all, for instance.
*/
//! return true if bodies bb overlap in all 3 dimensions
bool InsertionSortCollider::spatialOverlapPeri(Body::id_t id1, Body::id_t id2,Scene* scene, Vector2i& periods) const {
	assert(periodic);
	assert(id1!=id2); // programming error, or weird bodies (too large?)
	for(int axis=0; axis<2; axis++){
		Real dim=scene->cell->getSize()[axis];
		// LOG_DEBUG("dim["<<axis<<"]="<<dim);
		// too big bodies
		if (!allowBiggerThanPeriod){ assert(maxima[2*id1+axis]-minima[2*id1+axis]<.99*dim); assert(maxima[2*id2+axis]-minima[2*id2+axis]<.99*dim);}
    // define normalized positions relative to id1->max, and with +1 shift for id1->min so that body1's bounds cover an interval [shiftedMin; 1] at the end of a b1-centric period
		Real lmin = (minima[2*id2+axis]-maxima[2*id1+axis])*invSizes[axis];
		Real lmax = (maxima[2*id2+axis]-maxima[2*id1+axis])*invSizes[axis];
		Real shiftedMin = (minima[2*id1+axis]-maxima[2*id1+axis])*invSizes[axis]+1.;
		if((lmax-lmin)>0.99 || shiftedMin<0){
			if (allowBiggerThanPeriod) {periods[axis]=0; continue;}
			else {
        //LOG_FATAL("max="<<maxima[2*id2+axis]<<"min="<<minima[2*id2+axis]<<"id1"<<id1<<" id2="<<id2<<"lmin="<<lmin<<"lmax="<<lmax<<"shiftedMin="<<shiftedMin<<"period="<<periods[axis]);
				LOG_FATAL("Body #"<<((lmax-lmin)>0.99?id2:id1)<<" spans over 0.99 the cell size "<<dim<<" (axis="<<axis<<", see flag allowBiggerThanPeriod)");
				throw runtime_error(__FILE__ ": Body larger than half of the cell size encountered.");}
		}
		int period1 = floor(lmax);
		//overlap around zero, on the "+" side
		if ((lmin-period1) <= 1.e-7) {periods[axis]=-period1; continue;}
		 //overlap around 1, on the "-" side
		if ((lmax-period1+1e-7) >= shiftedMin) {periods[axis]=-period1-1; continue;}
		// none of the above, exit
    return false;
    /*

		// find body of which minimum when taken as period start will make the gap smaller
		Real m1=minima[2*id1+axis],m2=minima[2*id2+axis];
		Real wMn=(cellWrapRel(m1,m2,m2+dim)<cellWrapRel(m2,m1,m1+dim)) ? m2 : m1;

		int pmn1,pmx1,pmn2,pmx2;
		Real mn1=cellWrap(m1,wMn,wMn+dim,pmn1), mx1=cellWrap(maxima[2*id1+axis],wMn,wMn+dim,pmx1);
		Real mn2=cellWrap(m2,wMn,wMn+dim,pmn2), mx2=cellWrap(maxima[2*id2+axis],wMn,wMn+dim,pmx2);

		if((pmn1!=pmx1) || (pmn2!=pmx2)){
			if (allowBiggerThanPeriod) {
				// If both bodies are bigger, we place them in the (0,0,0) period
				if((pmn1!=pmx1) && (pmn2!=pmx2)) {periods[axis]=0;}
				// else we define period with the position of the small body (we assume the big one sits in period (0,0,0), keep that in mind if velGrad(.,axis) is not a null vector)
				else {
					//FIXME: not sure what to do here...
// 					periods[axis]=(pmn1==pmx1)? pmn1 : -pmn2;
					periods[axis]=0;
// 					return true;
				}
			} else {
				Real span=(pmn1!=pmx1?mx1-mn1:mx2-mn2); if(span<0) span=dim-span;
        //there is a bug
        sudodem::ObjectIO::save("tmp.xml.bz2","scene",scene);
				LOG_FATAL("Body #"<<(pmn1!=pmx1?id1:id2)<<" spans over half of the cell size "<<dim<<" (axis="<<axis<<", min="<<(pmn1!=pmx1?mn1:mn2)<<", max="<<(pmn1!=pmx1?mx1:mx2)<<", span="<<span<<", see flag allowBiggerThanPeriod)");
				throw runtime_error(__FILE__ ": Body larger than half of the cell size encountered.");}
		}
		else {
			periods[axis]=(int)(pmn1-pmn2);
			if(!(mn1<=mx2 && mx1 >= mn2)) return false;}
  */
	}
	return true;
}

boost::python::tuple InsertionSortCollider::dumpBounds(){
  boost::python::list bl[2]; // 2 bound lists, inserted into the tuple at the end
	for(int axis=0; axis<2; axis++){
		VecBounds& V=BB[axis];
		if(periodic){
			for(long i=0; i<V.size; i++){
				long ii=V.norm(i); // start from the period boundary
				bl[axis].append(boost::python::make_tuple(V[ii].coord,(V[ii].flags.isMin?-1:1)*V[ii].id,V[ii].period));
			}
		} else {
			for(long i=0; i<V.size; i++){
				bl[axis].append(boost::python::make_tuple(V[i].coord,(V[i].flags.isMin?-1:1)*V[i].id));
			}
		}
	}
	return boost::python::make_tuple(bl[0],bl[1]);
}
