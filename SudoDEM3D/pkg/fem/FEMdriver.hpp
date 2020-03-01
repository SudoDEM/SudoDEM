/*************************************************************************
 Copyright (C) 2017 by Sway Zhao		                                 *
*  zhswee@gmail.com      				                            	 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/common/FieldApplier.hpp>
#include<sudodem/core/Interaction.hpp>
//#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/common/Callbacks.hpp>
#include<sudodem/pkg/dem/GlobalStiffnessTimeStepper.hpp>

#include "TriElement.hpp"

#ifdef SUDODEM_OPENMP
	#include<omp.h>
#endif

class FEMdriver : public FieldApplier{
	inline void cundallDamp1st(Vector3r& force, const Vector3r& vel);
	inline void cundallDamp2nd(const Real& dt, const Vector3r& vel, Vector3r& accel);
	Quaternionr DotQ(const Vector3r& angVel, const Quaternionr& Q);

	// compute linear and angular acceleration, respecting State::blockedDOFs
	Vector3r computeAccel(const Vector3r& force, const Real& mass, int blockedDOFs);
	Vector3r computeAngAccel(const Vector3r& torque, const Vector3r& inertia, int blockedDOFs);
    //TriElement force
    void applyNodalForces();
    void driveNodes();
	#ifdef SUDODEM_OPENMP
	void ensureSync(); bool syncEnsured;
	#endif
	Matrix3r dVelGrad;

		#ifdef SUDODEM_OPENMP
			vector<Real> threadMaxVelocitySq;
		#endif
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(FEMdriver,GlobalEngine,"Engine integrating newtonian motion equations.",
		((Real,damping,0.2,,"damping coefficient for Cundall's non viscous damping (see [Chareyre2005]_) [-]"))
		//((Real,rot_damping,0.8,,"rotational damping coefficient for Cundall's non viscous damping (see [Chareyre2005]_) [-]"))
		((bool,bending,false,,"is TriElement bending?"))
        ((bool,rotIncr,false,,"is TriElement bending?"))
		((Real,thickness,0,,"Element thickness"))
        ((Real,young,0,,"Young's modulus"))
        ((Real,nu,0,,"Poisson's ratio"))
		((Real,bendThickness,0,,"Element bend thickness"))
        //((bool,quiet_system_flag,false,,"the flag for quiet system"))
		((Vector3r,gravity,Vector3r::Zero(),,"Gravitational acceleration (effectively replaces GravityEngine)."))
		((Real,maxVelocitySq,NaN,,"store square of max. velocity, for informative purposes; computed again at every step. |yupdate|"))
		((Vector3r,prevCellSize,Vector3r(NaN,NaN,NaN),Attr::hidden,"cell size from previous step, used to detect change and find max velocity"))
		((bool,warnNoForceReset,true,,"Warn when forces were not resetted in this step by :yref:`ForceResetter`; this mostly points to :yref:`ForceResetter` being forgotten incidentally and should be disabled only with a good reason."))
		((int,mask,-1,,"If mask defined and the bitwise AND between mask and body`s groupMask gives 0, the body will not move/rotate. Velocities and accelerations will be calculated not paying attention to this parameter."))
		,
		/*ctor*/
			
			#ifdef SUDODEM_OPENMP
				threadMaxVelocitySq.resize(omp_get_max_threads()); syncEnsured=false;
			#endif
		,/*py*/
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(FEMdriver);

