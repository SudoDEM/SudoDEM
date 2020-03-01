// 2008, 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#include<time.h>
#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>

class PeriodicEngine:  public GlobalEngine {
	public:
		static Real getClock(){ timeval tp; gettimeofday(&tp,NULL); return tp.tv_sec+tp.tv_usec/1e6; }
		virtual ~PeriodicEngine() {}; // vtable
		virtual bool isActivated(){
			const Real& virtNow=scene->time;
			Real realNow=getClock();
			const long& iterNow=scene->iter;
			if (iterNow<iterLast) nDone=0;//handle O.resetTime(), all counters will be initialized again
			if((nDo<0 || nDone<nDo) &&
				((virtPeriod>0 && virtNow-virtLast>=virtPeriod) ||
				 (realPeriod>0 && realNow-realLast>=realPeriod) ||
				 (iterPeriod>0 && iterNow-iterLast>=iterPeriod))){
				realLast=realNow; virtLast=virtNow; iterLast=iterNow; nDone++;
				return true;
			}
			if(nDone==0){
				realLast=realNow; virtLast=virtNow; iterLast=iterNow; nDone++;
				if(initRun) return true;
				return false;
			}
			return false;
		}
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(PeriodicEngine,GlobalEngine,
		"Run Engine::action with given fixed periodicity real time (=wall clock time, computation time), \
		virtual time (simulation time), iteration number), by setting any of those criteria \
		(virtPeriod, realPeriod, iterPeriod) to a positive value. They are all negative (inactive)\
		by default.\n\n\
		\
		The number of times this engine is activated can be limited by setting nDo>0. If the number of activations \
		will have been already reached, no action will be called even if an active period has elapsed. \n\n\
		\
		If initRun is set (false by default), the engine will run when called for the first time; otherwise it will only \
		start counting period (realLast etc interal variables) from that point, but without actually running, and will run \
		only once a period has elapsed since the initial run. \n\n\
		\
		This class should not be used directly; rather, derive your own engine which you want to be run periodically. \n\n\
		\
		Derived engines should override Engine::action(), which will be called periodically. If the derived Engine \
		overrides also Engine::isActivated, it should also take in account return value from PeriodicEngine::isActivated, \
		since otherwise the periodicity will not be functional. \n\n\
		\
		Example with :yref:`PyRunner`, which derives from PeriodicEngine; likely to be encountered in python scripts:: \n\n\
		\
			PyRunner(realPeriod=5,iterPeriod=10000,command='print O.iter')	\n\n\
		\
		will print iteration number every 10000 iterations or every 5 seconds of wall clock time, whiever comes first since it was \
		last run.",
		((Real,virtPeriod,((void)"deactivated",0),,"Periodicity criterion using virtual (simulation) time (deactivated if <= 0)"))
		((Real,realPeriod,((void)"deactivated",0),,"Periodicity criterion using real (wall clock, computation, human) time (deactivated if <=0)"))
		((long,iterPeriod,((void)"deactivated",0),,"Periodicity criterion using step number (deactivated if <= 0)"))
		((long,nDo,((void)"deactivated",-1),,"Limit number of executions by this number (deactivated if negative)"))
		((bool,initRun,false,,"Run the first time we are called as well."))
		((Real,virtLast,0,,"Tracks virtual time of last run |yupdate|."))
		((Real,realLast,0,,"Tracks real time of last run |yupdate|."))
		((long,iterLast,0,,"Tracks step number of last run |yupdate|."))
		((long,nDone,0,,"Track number of executions (cummulative) |yupdate|.")),
		/* this will be put inside the ctor */ realLast=getClock();
	);
};
REGISTER_SERIALIZABLE(PeriodicEngine);
