/*************************************************************************
*  Copyright (C) 2016 by Zhswee                                          *
*  zhswee@gmail.com                                                      *
*  Engine for expanding polySuperEllipsoids                              *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<sudodem/pkg/common/BoundaryController.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/dem/Superquadrics.hpp>
#include<boost/array.hpp>
#include<sudodem/core/PartialEngine.hpp>
#include<sudodem/pkg/common/Facet.hpp>
#include<fstream>
#include<string>
#ifdef SUDODEM_OPENMP
	#include<omp.h>
#endif

class State;

class Scene;
class State;


class Expander : public BoundaryController
{
	private :
		bool  firstRun;
		std::vector<double> expandAlphaList;
		Real curExpandSlice;//at which step of the expansion
		public :
		virtual ~Expander();
		virtual void action();
		void initial();
		void expand();
		//void getStressStrain();//
		SUDODEM_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(Expander,BoundaryController,
		"An engine to expand polySuperEllipsoids when preparing a specimen using the particle-expansion method\n"
		,
    ((double,initExpandAlpha,0.1,,"the init expansion alpha."))
    ((short int,expandSlice,100,,"steps needing to expand a particle to the target"))
    ((double,preExpandAlpha,0.1,Attr::readonly,"the previous expansion alpha"))
		((int, cycleInterval, 100,,"trigger a new expansion after cycleInterval steps of running"))
		,
		/* extra initializers */
		,
		//
		,
   		/* constructor */
   		firstRun =true;
		,
		)
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Expander);
