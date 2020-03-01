/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once
#include<sudodem/pkg/common/Facet.hpp>
#include<sudodem/core/Body.hpp>
#include<sudodem/core/BodyContainer.hpp>


class STLImporter {
    public:
	vector<shared_ptr<Body> > import(const char*);
	DECLARE_LOGGER;
};

