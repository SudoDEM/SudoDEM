/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Bronek Kozicki                                  *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#include<sudodem/lib/base/Math.hpp>
#include <dlfcn.h>

#include<sudodem/lib/base/Logging.hpp>
//#include<sudodem/lib/base/Math.hpp>

class DynLibManager
{
	private :
		std::map<const std::string, void *> handles;
		bool autoUnload;

	public :
		DynLibManager ();
		~DynLibManager ();
		void addBaseDirectory(const std::string& dir);

		bool load(const std::string& libName);

		bool unload (const std::string& libName);
		bool isLoaded (const std::string& libName);
		bool unloadAll ();
		void setAutoUnload ( bool enabled );

    std::string lastError();
		DECLARE_LOGGER;

	private :
		bool closeLib(const std::string libName);
		bool error();
    std::string lastError_;
};


