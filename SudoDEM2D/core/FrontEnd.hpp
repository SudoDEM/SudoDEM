/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include "Omega.hpp"

#include<sudodem/lib/factory/Factorable.hpp>

class FrontEnd : public Factorable
{
	public :
		FrontEnd () {};
		virtual ~FrontEnd () {};

		virtual int run(int , char * []) { return -1;};
		// called before actually invoking it
		virtual bool available(){return false;}

	REGISTER_CLASS_AND_BASE(FrontEnd,Factorable);
};
REGISTER_FACTORABLE(FrontEnd);


