/*************************************************************************
*  Copyright (C) 2006 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include <sudodem/lib/base/Logging.hpp>
#include "ThreadRunner.hpp"
#include "ThreadWorker.hpp"

#include <boost/thread/thread.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

CREATE_LOGGER(ThreadRunner);

void ThreadRunner::run()
{
	// this is the body of execution of separate thread
	boost::mutex::scoped_lock lock(m_runmutex);
	try{
		workerThrew=false;
		while(looping()) {
			call();
			if(m_thread_worker->shouldTerminate()){ stop(); return; }
		}
	} catch (std::exception& e){
		LOG_FATAL("Exception occured: "<<std::endl<<e.what());
		workerException=std::exception(e); workerThrew=true;
		stop(); return;
	}
}

void ThreadRunner::call()
{
	// this is the body of execution of separate thread
	//
	// FIXME - if several threads are blocked here and waiting, and the
	// destructor is called we get a crash. This happens if some other
	// thread is calling spwanSingleAction in a loop (instead of calling
	// start() and stop() as it normally should). This is currently the
	// case of SimulationController with synchronization turned on.
	//
	// the solution is to use a counter (perhaps recursive_mutex?) which
	// will count the number of threads in the queue, and only after they
	// all finish execution the destructor will be able to finish its work
	//
	boost::mutex::scoped_lock lock(m_callmutex);
	m_thread_worker->setTerminate(false);
	m_thread_worker->callSingleAction();
}

void ThreadRunner::pleaseTerminate()
{
	stop();
	m_thread_worker->setTerminate(true);
}

void ThreadRunner::spawnSingleAction()
{
	boost::mutex::scoped_lock boollock(m_boolmutex);
	boost::mutex::scoped_lock calllock(m_callmutex);
	if(m_looping) return;
	boost::function0<void> call( boost::bind( &ThreadRunner::call , this ) );
	boost::thread th(call);
}

void ThreadRunner::start()
{
	boost::mutex::scoped_lock lock(m_boolmutex);
	if(m_looping) return;
	m_looping=true;
	boost::function0<void> run( boost::bind( &ThreadRunner::run , this ) );
	boost::thread th(run);
}

void ThreadRunner::stop()
{
	//std::cerr<<__FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<std::endl;
	if(!m_looping) return;
	//std::cerr<<__FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<std::endl;
	boost::mutex::scoped_lock lock(m_boolmutex);
	//std::cerr<<__FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<std::endl;
	m_looping=false;
}

bool ThreadRunner::looping()
{
	boost::mutex::scoped_lock lock(m_boolmutex);
	return m_looping;
}

ThreadRunner::~ThreadRunner()
{
	pleaseTerminate();
	boost::mutex::scoped_lock runlock(m_runmutex);
	boost::mutex::scoped_lock calllock(m_callmutex);
}

