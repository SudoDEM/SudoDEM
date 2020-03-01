/*************************************************************************
*  Copyright (C) 2006 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include "Scene.hpp"
#include <sudodem/lib/base/Logging.hpp>
#include "SimulationFlow.hpp"
//#include "Scene.hpp"
#include "Omega.hpp"

CREATE_LOGGER(SimulationFlow);

void SimulationFlow::singleAction()
{
	Scene* scene=Omega::instance().getScene().get();
	if (!scene) throw logic_error("SimulationFlow::singleAction: no Scene object?!");
	if(scene->subStepping) { LOG_INFO("Sub-stepping disabled when running simulation continuously."); scene->subStepping=false; }
	scene->moveToNextTimeStep();
	if(scene->stopAtIter>0 && scene->iter==scene->stopAtIter) setTerminate(true);
	if(scene->stopAtTime>0 && scene->time==scene->stopAtTime) setTerminate(true);
};

