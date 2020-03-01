// Sway (C) 2017 zhswee@gmail.com
#pragma once

#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/fem/Node.hpp>

#include<boost/static_assert.hpp>
// make sure that (void*)&vec[0]==(void*)&vec
BOOST_STATIC_ASSERT(sizeof(Vector3r)==3*sizeof(Real));


#ifdef SUDODEM_OPENMP

#include<omp.h>
/*! Container for Node External Variables (forces), typically forces and torques from interactions.
 * Values should be reset at every iteration by calling NodeForceContainer::reset();
 * If you want to add your own force type, you need to:
 *
 * 	1. Create storage vector
 * 	2. Create accessor function
 * 	3. Update the resize function
 * 	4. Update the reset function
 * 	5. update the sync function (for the multithreaded implementation)
 *
 * This class exists in two flavors: non-parallel and parallel. The parallel one stores
 * force increments separately for every thread and sums those when sync() is called.
 * The reason of this design is that the container is not truly random-access, but rather
 * is written to everywhere in one phase and read in the next one. Adding to force/torque
 * marks the container as dirty and sync() must be performed before reading the stored data.
 * Calling getForce/getTorque when the container is not synchronized throws an exception.
 *
 * It is intentional that sync() needs to be called exlicitly, since syncs are expensive and
 * the programmer should be aware of that. Sync is however performed only if the container
 * is dirty. Every full sync increments the syncCount variable, that should ideally equal
 * the number of steps (one per step).
 *
 * The number of threads (omp_get_max_threads) may not change once NodeForceContainer is constructed.
 *
 * The non-parallel flavor has the same interface, but sync() is no-op and synchronization
 * is not enforced at all.
 */

//! This is the parallel flavor of NodeForceContainer
class NodeForceContainer{
	private:
		typedef std::vector<Vector3r> vvector;
		std::vector<vvector> _forceData;
		std::vector<vvector> _torqueData;
		std::vector<vvector> _moveData;
		std::vector<vvector> _rotData;
		std::vector<Node::id_t>  _maxId;
		vvector _force, _torque, _move, _rot, _permForce, _permTorque;
		std::vector<size_t> sizeOfThreads;
		size_t size;
		size_t permSize;
		bool syncedSizes;
		int nThreads;
		bool synced,moveRotUsed,permForceUsed;
		boost::mutex globalMutex;
		Vector3r _zero;

		inline void ensureSize(Node::id_t id, int threadN){
			assert(nThreads>omp_get_thread_num());
			const Node::id_t idMaxTmp = max(id, _maxId[threadN]);
			_maxId[threadN] = 0;
			if (threadN<0) {
				resizePerm(min((size_t)1.5*(idMaxTmp+100),(size_t)(idMaxTmp+2000)));
			} else if (sizeOfThreads[threadN]<=(size_t)idMaxTmp) {
				resize(min((size_t)1.5*(idMaxTmp+100),(size_t)(idMaxTmp+2000)),threadN);
			}
		}
		inline void ensureSynced(){ if(!synced) throw runtime_error("NodeForceContainer not thread-synchronized; call sync() first!"); }

		// dummy function to avoid template resolution failure
		friend class boost::serialization::access; template<class ArchiveT> void serialize(ArchiveT & ar, unsigned int version){}
	public:
		NodeForceContainer(): size(0), permSize(0),syncedSizes(true),synced(true),moveRotUsed(false),permForceUsed(false),_zero(Vector3r::Zero()),syncCount(0),lastReset(0){
			nThreads=omp_get_max_threads();
			for(int i=0; i<nThreads; i++){
				_forceData.push_back(vvector()); _torqueData.push_back(vvector());
				_moveData.push_back(vvector());  _rotData.push_back(vvector());
				sizeOfThreads.push_back(0);
				_maxId.push_back(0);
			}
		}
		const Vector3r& getForce(Node::id_t id)         { ensureSynced(); return ((size_t)id<size)?_force[id]:_zero; }
		void  addForce(Node::id_t id, const Vector3r& f){ ensureSize(id,omp_get_thread_num()); synced=false;   _forceData[omp_get_thread_num()][id]+=f;}
		const Vector3r& getTorque(Node::id_t id)        { ensureSynced(); return ((size_t)id<size)?_torque[id]:_zero; }
		void addTorque(Node::id_t id, const Vector3r& t){ ensureSize(id,omp_get_thread_num()); synced=false;   _torqueData[omp_get_thread_num()][id]+=t;}
		const Vector3r& getMove(Node::id_t id)          { ensureSynced(); return ((size_t)id<size)?_move[id]:_zero; }
		void  addMove(Node::id_t id, const Vector3r& m) { ensureSize(id,omp_get_thread_num()); synced=false; moveRotUsed=true; _moveData[omp_get_thread_num()][id]+=m;}
		const Vector3r& getRot(Node::id_t id)           { ensureSynced(); return ((size_t)id<size)?_rot[id]:_zero; }
		void  addRot(Node::id_t id, const Vector3r& r)  { ensureSize(id,omp_get_thread_num()); synced=false; moveRotUsed=true; _rotData[omp_get_thread_num()][id]+=r;}
		void  addMaxId(Node::id_t id)                   { _maxId[omp_get_thread_num()]=id;}

		void  addPermForce(Node::id_t id, const Vector3r& f){ ensureSize(id,-1); synced=false;   _permForce[id]=f; permForceUsed=true;}
		void addPermTorque(Node::id_t id, const Vector3r& t){ ensureSize(id,-1); synced=false;   _permTorque[id]=t; permForceUsed=true;}
		const Vector3r& getPermForce(Node::id_t id) { ensureSynced(); return ((size_t)id<size)?_permForce[id]:_zero; }
		const Vector3r& getPermTorque(Node::id_t id) { ensureSynced(); return ((size_t)id<size)?_permTorque[id]:_zero; }

		/*! Function to allow friend classes to get force even if not synced. Used for clumps by NewtonIntegrator.
		* Dangerous! The caller must know what it is doing! (i.e. don't read after write
		* for a particular Node id. */
		const Vector3r& getForceUnsynced (Node::id_t id){assert ((size_t)id<size); return _force[id];}
		const Vector3r& getTorqueUnsynced(Node::id_t id){assert ((size_t)id<size); return _torque[id];}
		void  addForceUnsynced(Node::id_t id, const Vector3r& f){ assert ((size_t)id<size); _force[id]+=f; }
		void  addTorqueUnsynced(Node::id_t id, const Vector3r& m){ assert ((size_t)id<size); _torque[id]+=m; }

		/* To be benchmarked: sum thread data in getForce/getTorque upon request for each Node individually instead of by the sync() function globally */
		// this function is used from python so that running simulation is not slowed down by sync'ing on occasions
		// since Vector3r writes are not atomic, it might (rarely) return wrong value, if the computation is running meanwhile
		Vector3r getForceSingle (Node::id_t id){ Vector3r ret(Vector3r::Zero()); for(int t=0; t<nThreads; t++){ ret+=((size_t)id<sizeOfThreads[t])?_forceData [t][id]:_zero; } if (permForceUsed) ret+=_permForce[id]; return ret; }
		Vector3r getTorqueSingle(Node::id_t id){ Vector3r ret(Vector3r::Zero()); for(int t=0; t<nThreads; t++){ ret+=((size_t)id<sizeOfThreads[t])?_torqueData[t][id]:_zero; } if (permForceUsed) ret+=_permTorque[id]; return ret; }
		Vector3r getMoveSingle  (Node::id_t id){ Vector3r ret(Vector3r::Zero()); for(int t=0; t<nThreads; t++){ ret+=((size_t)id<sizeOfThreads[t])?_moveData  [t][id]:_zero; } return ret; }
		Vector3r getRotSingle   (Node::id_t id){ Vector3r ret(Vector3r::Zero()); for(int t=0; t<nThreads; t++){ ret+=((size_t)id<sizeOfThreads[t])?_rotData   [t][id]:_zero; } return ret; }

		inline void syncSizesOfContainers() {
			if (syncedSizes) return;
			//check whether all containers have equal length, and if not resize it
			for(int i=0; i<nThreads; i++){
				if (sizeOfThreads[i]<size) resize(size,i);
			}
			_force.resize(size,Vector3r::Zero());
			_torque.resize(size,Vector3r::Zero());
			_permForce.resize(size,Vector3r::Zero());
			_permTorque.resize(size,Vector3r::Zero());
			_move.resize(size,Vector3r::Zero());
			_rot.resize(size,Vector3r::Zero());
			syncedSizes=true;
		}
		/* Sum contributions from all threads, save to _force&_torque.
		 * Locks globalMutex, since one thread modifies common data (_force&_torque).
		 * Must be called before get* methods are used. Exception is thrown otherwise, since data are not consistent. */
		inline void sync(){
			for(int i=0; i<nThreads; i++){
				if (_maxId[i] > 0) { synced = false;}
			}
			if(synced) return;
			boost::mutex::scoped_lock lock(globalMutex);
			if(synced) return; // if synced meanwhile

			for(int i=0; i<nThreads; i++){
				if (_maxId[i] > 0) { ensureSize(_maxId[i],i);}
			}

			syncSizesOfContainers();

			for(long id=0; id<(long)size; id++){
				Vector3r sumF(Vector3r::Zero()), sumT(Vector3r::Zero());
				for(int thread=0; thread<nThreads; thread++){ sumF+=_forceData[thread][id]; sumT+=_torqueData[thread][id];}
				_force[id]=sumF; _torque[id]=sumT;
				if (permForceUsed) {_force[id]+=_permForce[id]; _torque[id]+=_permTorque[id];}
			}
			if(moveRotUsed){
				for(long id=0; id<(long)size; id++){
					Vector3r sumM(Vector3r::Zero()), sumR(Vector3r::Zero());
					for(int thread=0; thread<nThreads; thread++){ sumM+=_moveData[thread][id]; sumR+=_rotData[thread][id];}
					_move[id]=sumM; _rot[id]=sumR;
				}
			}
			synced=true; syncCount++;
		}
		unsigned long syncCount;
		long lastReset;

		void resize(size_t newSize, int threadN){
			_forceData [threadN].resize(newSize,Vector3r::Zero());
			_torqueData[threadN].resize(newSize,Vector3r::Zero());
			_moveData[threadN].resize(newSize,Vector3r::Zero());
			_rotData[threadN].resize(newSize,Vector3r::Zero());
			sizeOfThreads[threadN] = newSize;
			if (size<newSize) size=newSize;
			syncedSizes=false;
		}
		void resizePerm(size_t newSize){
			_permForce.resize(newSize,Vector3r::Zero());
			_permTorque.resize(newSize,Vector3r::Zero());
			permSize = newSize;
			if (size<newSize) size=newSize;
			syncedSizes=false;
		}
		/*! Reset all resetable data, also reset summary forces/torques and mark the container clean.
		If resetAll, reset also user defined forces and torques*/
		// perhaps should be private and friend Scene or whatever the only caller should be
		void reset(long iter, bool resetAll=false){
			syncSizesOfContainers();
			for(int thread=0; thread<nThreads; thread++){
				memset(&_forceData [thread][0],0,sizeof(Vector3r)*sizeOfThreads[thread]);
				memset(&_torqueData[thread][0],0,sizeof(Vector3r)*sizeOfThreads[thread]);
				if(moveRotUsed){
					memset(&_moveData  [thread][0],0,sizeof(Vector3r)*sizeOfThreads[thread]);
					memset(&_rotData   [thread][0],0,sizeof(Vector3r)*sizeOfThreads[thread]);
				}
			}
			memset(&_force [0], 0,sizeof(Vector3r)*size);
			memset(&_torque[0], 0,sizeof(Vector3r)*size);
			if(moveRotUsed){
				memset(&_move  [0], 0,sizeof(Vector3r)*size);
				memset(&_rot   [0], 0,sizeof(Vector3r)*size);
			}
			if (resetAll){
				memset(&_permForce [0], 0,sizeof(Vector3r)*size);
				memset(&_permTorque[0], 0,sizeof(Vector3r)*size);
				permForceUsed = false;
			}
			if (!permForceUsed) synced=true; else synced=false;
			moveRotUsed=false;
			lastReset=iter;
		}
		//! say for how many threads we have allocated space
		const int& getNumAllocatedThreads() const {return nThreads;}
		const bool& getMoveRotUsed() const {return moveRotUsed;}
		const bool& getPermForceUsed() const {return permForceUsed;}
};

#else
//! This is the non-parallel flavor of NodeForceContainer
class NodeForceContainer {
	private:
		std::vector<Vector3r> _force;
		std::vector<Vector3r> _torque;
		std::vector<Vector3r> _move;
		std::vector<Vector3r> _rot;
		std::vector<Vector3r> _permForce, _permTorque;
		Node::id_t _maxId;
		size_t size;
		size_t permSize;
		inline void ensureSize(Node::id_t id){
			const Node::id_t idMaxTmp = max(id, _maxId);
			_maxId = 0;
			if(size<=(size_t)idMaxTmp) resize(min((size_t)1.5*(idMaxTmp+100),(size_t)(idMaxTmp+2000)));
		}
		#if 0
			const Vector3r& getForceUnsynced (Node::id_t id){ return getForce(id);}
			const Vector3r& getTorqueUnsynced(Node::id_t id){ return getForce(id);}
		#endif
		bool moveRotUsed, permForceUsed;
		// dummy function to avoid template resolution failure
		friend class boost::serialization::access; template<class ArchiveT> void serialize(ArchiveT & ar, unsigned int version){}
	public:
		NodeForceContainer(): size(0), permSize(0), moveRotUsed(false), permForceUsed(false), syncCount(0), lastReset(0), _maxId(0){}
		const Vector3r& getForce(Node::id_t id){ensureSize(id); return _force[id];}
		void  addForce(Node::id_t id,const Vector3r& f){ensureSize(id); _force[id]+=f;}
		const Vector3r& getTorque(Node::id_t id){ensureSize(id); return _torque[id];}
		void  addTorque(Node::id_t id,const Vector3r& t){ensureSize(id); _torque[id]+=t;}
		const Vector3r& getMove(Node::id_t id){ensureSize(id); return _move[id];}
		void  addMove(Node::id_t id,const Vector3r& f){ensureSize(id); moveRotUsed=true; _move[id]+=f;}
		const Vector3r& getRot(Node::id_t id){ensureSize(id); return _rot[id];}
		void  addRot(Node::id_t id,const Vector3r& f){ensureSize(id); moveRotUsed=true; _rot[id]+=f;}
		void  addPermForce(Node::id_t id, const Vector3r& f){ ensureSize(id);  _permForce[id]=f; permForceUsed=true;}
		void addPermTorque(Node::id_t id, const Vector3r& t){ ensureSize(id);  _permTorque[id]=t; permForceUsed=true;}
		void  addMaxId(Node::id_t id) { _maxId=id;}
		const Vector3r& getPermForce(Node::id_t id) { ensureSize(id); return _permForce[id]; }
		const Vector3r& getPermTorque(Node::id_t id) { ensureSize(id); return _permTorque[id]; }
		// single getters do the same as globally synced ones in the non-parallel flavor
		const Vector3r getForceSingle (Node::id_t id){
			ensureSize(id);
			if (permForceUsed) {
				return _force [id] + _permForce[id];
			} else {
				return _force [id];
			}
		}
		const Vector3r getTorqueSingle(Node::id_t id){
			ensureSize(id);
			if (permForceUsed) {
				return _torque[id] + _permTorque[id];
			} else {
				return _torque[id];
			}
		}
		const Vector3r& getMoveSingle  (Node::id_t id){ ensureSize(id); return _move  [id]; }
		const Vector3r& getRotSingle   (Node::id_t id){ ensureSize(id); return _rot   [id]; }

		//! Set all forces to zero
		void reset(long iter, bool resetAll=false){
			memset(&_force [0],0,sizeof(Vector3r)*size);
			memset(&_torque[0],0,sizeof(Vector3r)*size);
			if(moveRotUsed){
				memset(&_move  [0],0,sizeof(Vector3r)*size);
				memset(&_rot   [0],0,sizeof(Vector3r)*size);
				moveRotUsed=false;
			}
			if (resetAll){
				memset(&_permForce [0], 0,sizeof(Vector3r)*size);
				memset(&_permTorque[0], 0,sizeof(Vector3r)*size);
				permForceUsed = false;
			}
			lastReset=iter;
		}

		void sync() {
			if (_maxId>0) {ensureSize(_maxId); _maxId=0;}
			if (permForceUsed) {
				for(long id=0; id<(long)size; id++) {
					_force[id]+=_permForce[id];
					_torque[id]+=_permTorque[id];
				}
			}
			return;
		}
		unsigned long syncCount;
		// interaction in which the container was last reset; used by NewtonIntegrator to detect whether ForceResetter was not forgotten
		long lastReset;
		/*! Resize the container; this happens automatically,
		 * but you may want to set the size beforehand to avoid resizes as the simulation grows. */
		void resize(size_t newSize){
			_force.resize(newSize,Vector3r::Zero());
			_torque.resize(newSize,Vector3r::Zero());
			_permForce.resize(newSize,Vector3r::Zero());
			_permTorque.resize(newSize,Vector3r::Zero());
			_move.resize(newSize,Vector3r::Zero());
			_rot.resize(newSize,Vector3r::Zero());
			size=newSize;
		}
		const int getNumAllocatedThreads() const {return 1;}
		const bool& getMoveRotUsed() const {return moveRotUsed;}
		const bool& getPermForceUsed() const {return permForceUsed;}
};

#endif
