/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  with help from Bronek Kozicki                                         *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "DynLibManager.hpp"
#include "ClassFactory.hpp"

CREATE_LOGGER(DynLibManager);


DynLibManager::DynLibManager ()
{
	autoUnload = true;
}


DynLibManager::~DynLibManager ()
{
	if(autoUnload) unloadAll();
}

// load plugin with given filename
bool DynLibManager::load (const string& lib){
	if (lib.empty()) throw std::runtime_error(__FILE__ ": got empty library name to load.");
	void* handle = dlopen(lib.c_str(),RTLD_GLOBAL | RTLD_NOW);
	if (!handle) return !error();
	handles[lib] = handle;
	return true;
}

// unload plugin, given full filename
bool DynLibManager::unload (const string& libName)
{
	if (isLoaded(libName))
	return closeLib(libName);
	else return false;
}


bool DynLibManager::unloadAll ()
{
	std::map<const string, void *>::iterator ith  = handles.begin();
	std::map<const string, void *>::iterator ithEnd  = handles.end();
	for( ; ith!=ithEnd ; ++ith)
		if ((*ith).first.length()!=0)
			unload((*ith).first);
	return false;
}


bool DynLibManager::isLoaded (const string& libName)
{
	std::map<const string, void *>::iterator ith = handles.find(libName);	
	return (ith!= handles.end() && (*ith).second!=NULL);
}


void DynLibManager::setAutoUnload ( bool enabled )
{
	autoUnload = enabled;
}


bool DynLibManager::closeLib(const string libName)
{
	dlclose(handles[libName]);
	return !error();

}

std::string DynLibManager::lastError()
{
        return lastError_;
}

bool DynLibManager::error() 
{
 	char * error = dlerror();
	if (error != NULL)  { lastError_ = error;	}
	return (error!=NULL);
}


