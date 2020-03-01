/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<sudodem/lib/serialization/Serializable.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Timing.hpp>
#include<sudodem/lib/base/Logging.hpp>

class Body;
class Scene;

CREATE_LOGGER(Engine);

class Engine: public Serializable{
	public:
		// pointer to the simulation, set at every step by Scene::moveToNextTimeStep
		Scene* scene;
		//! high-level profiling information; not serializable
		TimingInfo timingInfo;
		//! precise profiling information (timing of fragments of the engine)
		shared_ptr<TimingDeltas> timingDeltas;
		virtual ~Engine() {};

		virtual bool isActivated() { return true; };
		virtual void action() {
			LOG_FATAL("Engine "<<getClassName()<<" calling virtual method Engine::action(). Please submit bug report at http://bugs.launchpad.net/sudodem.");
			throw std::logic_error("Engine::action() called.");
		}
	private:
		// py access funcs
		TimingInfo::delta timingInfo_nsec_get(){return timingInfo.nsec;};
		void timingInfo_nsec_set(TimingInfo::delta d){ timingInfo.nsec=d;}
		long timingInfo_nExec_get(){return timingInfo.nExec;};
		void timingInfo_nExec_set(long d){ timingInfo.nExec=d;}
		void explicitAction() {scene=Omega::instance().getScene().get(); action();};

	DECLARE_LOGGER;

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Engine,Serializable,"Basic execution unit of simulation, called from the simulation loop (O.engines)",
		((bool,dead,false,,"If true, this engine will not run at all; can be used for making an engine temporarily deactivated and only resurrect it at a later point."))
		((int, ompThreads, -1,,"Number of threads to be used in the engine. If ompThreads<0 (default), the number will be typically OMP_NUM_THREADS or the number N defined by 'sudodem -jN' (this behavior can depend on the engine though). This attribute will only affect engines whose code includes openMP parallel regions (e.g. :yref:`InteractionLoop`). This attribute is mostly useful for experiments or when combining :yref:`ParallelEngine` with engines that run parallel regions, resulting in nested OMP loops with different number of threads at each level."))
		((string,label,,,"Textual label for this object; must be valid python identifier, you can refer to it directly from python.")),
		/* ctor */ scene=Omega::instance().getScene().get();
		#ifdef USE_TIMING_DELTAS
			timingDeltas=shared_ptr<TimingDeltas>(new TimingDeltas);
		#endif
		,
		/* py */
		.add_property("execTime",&Engine::timingInfo_nsec_get,&Engine::timingInfo_nsec_set,"Cummulative time this Engine took to run (only used if :yref:`O.timingEnabled<Omega.timingEnabled>`\\ ==\\ ``True``).")
		.add_property("execCount",&Engine::timingInfo_nExec_get,&Engine::timingInfo_nExec_set,"Cummulative count this engine was run (only used if :yref:`O.timingEnabled<Omega.timingEnabled>`\\ ==\\ ``True``).")
		.def_readonly("timingDeltas",&Engine::timingDeltas,"Detailed information about timing inside the Engine itself. Empty unless enabled in the source code and :yref:`O.timingEnabled<Omega.timingEnabled>`\\ ==\\ ``True``.")
		.def("__call__",&Engine::explicitAction)
	);
};
REGISTER_SERIALIZABLE(Engine);



