/*************************************************************************
*  Copyright (C) 2006 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <boost/thread/mutex.hpp>
#include <boost/any.hpp>

class ThreadRunner;

/*! 
\brief	ThreadWorker contains information about tasks to be performed when
	the separate thread is executed.
 */

class ThreadWorker	//         perhaps simulation steps, or stage? as it is a single stage
			//         of the simulation, that consists of several steps
			// Update: it is more general now. simulation stages perhaps will be derived from this class
{
	private:
		/// You should check out ThreadRunner, it is used for execution control of this class
		friend class ThreadRunner;
		bool		m_should_terminate;
		bool		m_done;
		boost::mutex	m_mutex;
		boost::any	m_val;
		float		m_progress;
		std::string	m_status;
		void		callSingleAction();

	protected:
		void		setTerminate(bool);
		/// singleAction() can check whether someone asked for termination, and terminate if/when possible
		bool		shouldTerminate();
		/// if something must be returned, set the result using this method
		void		setReturnValue(boost::any);
		/// if you feel monitored for progress, you can set it here: a value between 0.0 and 1.0
		void		setProgress(float);
		/// if you feel being monitored for what currently is done, set the message here
		void		setStatus(std::string);
		/// derived classes must define this method, that's what is executed in separate thread
		virtual void	singleAction() = 0;

	public:
		ThreadWorker() : m_should_terminate(false), m_done(false), m_progress(0) {};
		virtual		~ThreadWorker() {};

		/// Returns a value between 0.0 and 1.0. Useful for updating a progress bar.
		float		progress(); // get_progress ? (pick a naming convention, efngh)
		/// You can display a message in GUI about what is the current work status
		std::string	getStatus();
		/// Check whether execution is finished,
		bool		done();
		/// then get the result.
		boost::any	getReturnValue();
};


