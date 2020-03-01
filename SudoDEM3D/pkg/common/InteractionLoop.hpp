// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/pkg/common/Callbacks.hpp>
#include<sudodem/pkg/common/Dispatching.hpp>

#ifdef USE_TIMING_DELTAS
	#define TIMING_DELTAS_CHECKPOINT(cpt) timingDeltas->checkpoint(cpt)
	#define TIMING_DELTAS_START() timingDeltas->start()
#else
	#define TIMING_DELTAS_CHECKPOINT(cpt)
	#define TIMING_DELTAS_START()
#endif

class InteractionLoop: public GlobalEngine {
	bool alreadyWarnedNoCollider;
	typedef std::pair<Body::id_t, Body::id_t> idPair;
	// store interactions that should be deleted after loop in action, not later
	#ifdef SUDODEM_OPENMP
		std::vector<std::list<idPair> > eraseAfterLoopIds;
		void eraseAfterLoop(Body::id_t id1,Body::id_t id2){ eraseAfterLoopIds[omp_get_thread_num()].push_back(idPair(id1,id2)); }
	#else
		list<idPair> eraseAfterLoopIds;
		void eraseAfterLoop(Body::id_t id1,Body::id_t id2){ eraseAfterLoopIds.push_back(idPair(id1,id2)); }
	#endif
	public:
		virtual void pyHandleCustomCtorArgs(boost::python::tuple& t, boost::python::dict& d);
		virtual void action();
		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(InteractionLoop,GlobalEngine,"Unified dispatcher for handling interaction loop at every step, for parallel performance reasons.\n\n.. admonition:: Special constructor\n\n\tConstructs from 3 lists of :yref:`Ig2<IGeomFunctor>`, :yref:`Ip2<IPhysFunctor>`, :yref:`Law<LawFunctor>` functors respectively; they will be passed to interal dispatchers, which you might retrieve.",
			((shared_ptr<IGeomDispatcher>,geomDispatcher,new IGeomDispatcher,Attr::readonly,":yref:`IGeomDispatcher` object that is used for dispatch."))
			((shared_ptr<IPhysDispatcher>,physDispatcher,new IPhysDispatcher,Attr::readonly,":yref:`IPhysDispatcher` object used for dispatch."))
			((shared_ptr<LawDispatcher>,lawDispatcher,new LawDispatcher,Attr::readonly,":yref:`LawDispatcher` object used for dispatch."))
			((vector<shared_ptr<IntrCallback> >,callbacks,,,":yref:`Callbacks<IntrCallback>` which will be called for every :yref:`Interaction`, if activated."))
			((bool, eraseIntsInLoop, false,,"Defines if the interaction loop should erase pending interactions, else the collider takes care of that alone (depends on what collider is used)."))
			,
			/*ctor*/ alreadyWarnedNoCollider=false;
				#ifdef SUDODEM_OPENMP
					eraseAfterLoopIds.resize(omp_get_max_threads());
				#endif
			,
			/*py*/
		);
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(InteractionLoop);
