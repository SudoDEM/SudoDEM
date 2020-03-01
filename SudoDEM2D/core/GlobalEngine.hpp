/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include "Engine.hpp"

class GlobalEngine: public Engine{
	public :
		virtual ~GlobalEngine() {};
	SUDODEM_CLASS_BASE_DOC(GlobalEngine,Engine,"Engine that will generally affect the whole simulation (contrary to PartialEngine).");
};
REGISTER_SERIALIZABLE(GlobalEngine);


