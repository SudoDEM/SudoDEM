/*************************************************************************
*  Copyright (C) 2008 by Jerome Duriez                                   *
*  jerome.duriez@hmg.inpg.fr                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<sudodem/core/Omega.hpp>
#include<sudodem/pkg/common/BoundaryController.hpp>
#include<sudodem/core/Body.hpp>



class Disp2DPropLoadEngine : public BoundaryController
{
	private :
		Real	dgamma	// the increment of horizontal displacement in one timestep, part of the disturbation
			,dh	// the increment of vertical displacement in one timestep, part of the disturbation
			,H0	// the height of the top box, at the beginnig of the application of the disturbation
			,X0	// the X-position of the top box, at the beginnig of the application of the disturbation
			,Fn0,Ft0// the normal and tangential force acting on the top box, at...
			,coordSs0,coordTot0
			;
		std::ofstream ofile;

		Real	alpha	// angle from the lower plate to the left box (trigo wise), as in other shear Engines, but here the Engine is able to find itself its value !
			,dalpha	// the increment over alpha
			;

		bool	firstIt;// true if this is the first iteration, false else.

		int	it_begin// the number of the it at which the computation starts
			;

		shared_ptr<Body> leftbox;
		shared_ptr<Body> rightbox;
		shared_ptr<Body> frontbox;
		shared_ptr<Body> backbox;
		shared_ptr<Body> topbox;
		shared_ptr<Body> boxbas;
		void saveData();
		void letDisturb();
		void stopMovement();		// to cancel all the velocities


	public :
		void 	action()
			,computeAlpha()
			;
		void postLoad(Disp2DPropLoadEngine&);

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(Disp2DPropLoadEngine,BoundaryController,
		"Disturbs a simple shear sample in a given displacement direction\n\nThis engine allows one to apply, on a simple shear sample, a loading controlled by du/dgamma = cste, which is equivalent to du + cste' * dgamma = 0 (proportionnal path loadings).\nTo do so, the upper plate of the simple shear box is moved in a given direction (corresponding to a given du/dgamma), whereas lateral plates are moved so that the box remains closed.\nThis engine can easily be used to perform directionnal probes, with a python script launching successivly the same .xml which contains this engine, after having modified the direction of loading (see *theta* attribute). That's why this Engine contains a *saveData* procedure which can save data on the state of the sample at the end of the loading (in case of successive loadings - for successive directions - through a python script, each line would correspond to one direction of loading).",
		((Body::id_t,id_topbox,3,,"the id of the upper wall"))
		((Body::id_t,id_boxbas,1,,"the id of the lower wall"))
		((Body::id_t,id_boxleft,0,,"the id of the left wall"))
		((Body::id_t,id_boxright,2,,"the id of the right wall"))
		((Body::id_t,id_boxfront,5,,"the id of the wall in front of the sample"))
		((Body::id_t,id_boxback,4,,"the id of the wall at the back of the sample"))
		((Real,theta,0.0,,"the angle, in a (gamma,h=-u) plane from the gamma - axis to the perturbation vector (trigo wise) [degrees]"))
		((Real,v,0.0,,"the speed at which the perturbation is imposed. In case of samples which are more sensitive to normal loadings than tangential ones, one possibility is to take v = V_shear - | (V_shear-V_comp)*sin(theta) | => v=V_shear in shear; V_comp in compression [m/s]"))
		((int,nbre_iter,0,,"the number of iterations of loading to perform"))
		((string,Key,"",,"string to add at the names of the saved files, and of the output file filled by *saveData*"))
		((bool,LOG,false,,"boolean controling the output of messages on the screen")),
		firstIt=true;
		alpha=Mathr::PI/2.0;
	)




};

REGISTER_SERIALIZABLE(Disp2DPropLoadEngine);

