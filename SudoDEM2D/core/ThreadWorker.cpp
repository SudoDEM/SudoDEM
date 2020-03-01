/*************************************************************************
*  Copyright (C) 2006 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "ThreadWorker.hpp"

void ThreadWorker::setTerminate(bool b)
{
	boost::mutex::scoped_lock lock(m_mutex);
	m_should_terminate=b;
};

bool ThreadWorker::shouldTerminate()
{
	boost::mutex::scoped_lock lock(m_mutex);
	return m_should_terminate;
};
		
void ThreadWorker::setProgress(float i)
{
	boost::mutex::scoped_lock lock(m_mutex);
	m_progress=i;
};

void ThreadWorker::setStatus(std::string s)
{
	boost::mutex::scoped_lock lock(m_mutex);
	m_status=s;
};

float ThreadWorker::progress()
{
	boost::mutex::scoped_lock lock(m_mutex);
	return m_progress;
};

std::string ThreadWorker::getStatus()
{
	boost::mutex::scoped_lock lock(m_mutex);
	return m_status;
};

void ThreadWorker::setReturnValue(boost::any a)
{
	boost::mutex::scoped_lock lock(m_mutex);
	m_val = a;
};

boost::any ThreadWorker::getReturnValue()
{
	boost::mutex::scoped_lock lock(m_mutex);
	return m_val;
};

bool ThreadWorker::done()
{
	boost::mutex::scoped_lock lock(m_mutex);
	return m_done;
};

void ThreadWorker::callSingleAction()
{
	{
		boost::mutex::scoped_lock lock(m_mutex);
		m_done = false;
	}
	this->singleAction();
	{
		boost::mutex::scoped_lock lock(m_mutex);
		m_done = true;
	}
};

