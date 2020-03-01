// 2008 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/PeriodicEngines.hpp>
#include<sudodem/lib/pyutil/gil.hpp>

class PyRunner: public PeriodicEngine {
	public :
		/* virtual bool isActivated: not overridden, PeriodicEngine handles that */
		virtual void action(){ if(command.size()>0) pyRunString(command); }
	SUDODEM_CLASS_BASE_DOC_ATTRS(PyRunner,PeriodicEngine,
		"Execute a python command periodically, with defined (and adjustable) periodicity. See :yref:`PeriodicEngine` documentation for details.",
		((string,command,"",,"Command to be run by python interpreter. Not run if empty."))
	);
};
REGISTER_SERIALIZABLE(PyRunner);
