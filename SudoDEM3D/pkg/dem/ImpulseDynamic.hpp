/*************************************************************************
 Copyright (C) 2018 by Shiwei Zhao		                         *
*  zhswee@gmail.com      					 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once
#ifdef __IBD__

#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/common/FieldApplier.hpp>
#include<sudodem/core/Interaction.hpp>
//#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/common/Callbacks.hpp>
#include<sudodem/pkg/dem/GlobalStiffnessTimeStepper.hpp>
#include"Superquadrics.hpp"

#ifdef SUDODEM_OPENMP
	#include<omp.h>
#endif

/*! An engine that can replace the usual series of engines used for integrating the laws of motion.

 */
class State;

class ImpulseDynamic : public FieldApplier{
	inline void cundallDamp1st(Vector3r& force, const Vector3r& vel);
	inline void cundallDamp2nd(const Real& dt, const Vector3r& vel, Vector3r& accel);
	inline void leapfrogTranslate(State*, const Body::id_t& id, const Real& dt); // leap-frog translate
	inline void leapfrogSphericalRotate(State*, const Body::id_t& id, const Real& dt); // leap-frog rotate of spherical body
	inline void leapfrogAsphericalRotate(State*, const Body::id_t& id, const Real& dt, const Vector3r& M); // leap-frog rotate of aspherical body
	inline void leapfrogSuperquadricsRotate(Superquadrics* ,State*, const Body::id_t& id, const Real& dt, const Vector3r& M); // leap-frog rotate of Superquadrics
	Quaternionr DotQ(const Vector3r& angVel, const Quaternionr& Q);

	// compute linear and angular acceleration, respecting State::blockedDOFs
	Vector3r computeAccel(const Vector3r& force, const Real& mass, int blockedDOFs);
	Vector3r computeAngAccel(const Vector3r& torque, const Vector3r& inertia, int blockedDOFs);

	void updateEnergy(const shared_ptr<Body>&b, const State* state, const Vector3r& fluctVel, const Vector3r& f, const Vector3r& m);
	#ifdef SUDODEM_OPENMP
	void ensureSync(); bool syncEnsured;
	#endif
	// whether the cell has changed from the previous step
	bool cellChanged;
	bool homoDeform;

	// wether a body has been selected in Qt view
	bool bodySelected;
	Matrix3r dVelGrad;

	public:
		bool densityScaling;// internal for density scaling
		Real updatingDispFactor;//(experimental) Displacement factor used to trigger bound update: the bound is updated only if updatingDispFactor*disp>sweepDist when >0, else all bounds are updated.
		// function to save maximum velocity, for the verlet-distance optimization
		void saveMaximaVelocity(const Body::id_t& id, State* state);
		void saveMaximaDisplacement(const shared_ptr<Body>& b);
		bool get_densityScaling ();
		void set_densityScaling (bool dsc);

		#ifdef SUDODEM_OPENMP
			vector<Real> threadMaxVelocitySq;
		#endif
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(ImpulseDynamic,GlobalEngine,"Engine integrating newtonian motion equations.",
		((Real,damping,0.2,,"damping coefficient for Cundall's non viscous damping (see [Chareyre2005]_) [-]"))
		//((Real,rot_damping,0.8,,"rotational damping coefficient for Cundall's non viscous damping (see [Chareyre2005]_) [-]"))
		((bool,isSuperquadrics,false,,"Enable optimation for Superquadrics"))//zhswee
        ((bool,quiet_system_flag,false,,"the flag for quiet system"))//zhswee
		((Vector3r,gravity,Vector3r::Zero(),,"Gravitational acceleration (effectively replaces GravityEngine)."))
		((Real,maxVelocitySq,NaN,,"store square of max. velocity, for informative purposes; computed again at every step. |yupdate|"))
		((bool,exactAsphericalRot,true,,"Enable more exact body rotation integrator for :yref:`aspherical bodies<Body.aspherical>` *only*, using formulation from [Allen1989]_, pg. 89."))
		((Matrix3r,prevVelGrad,Matrix3r::Zero(),,"Store previous velocity gradient (:yref:`Cell::velGrad`) to track acceleration. |yupdate|"))
		#ifdef SUDODEM_BODY_CALLBACK
			((vector<shared_ptr<BodyCallback> >,callbacks,,,"List (std::vector in c++) of :yref:`BodyCallbacks<BodyCallback>` which will be called for each body as it is being processed."))
		#endif
		((Vector3r,prevCellSize,Vector3r(NaN,NaN,NaN),Attr::hidden,"cell size from previous step, used to detect change and find max velocity"))
		((bool,warnNoForceReset,true,,"Warn when forces were not resetted in this step by :yref:`ForceResetter`; this mostly points to :yref:`ForceResetter` being forgotten incidentally and should be disabled only with a good reason."))
		// energy tracking
		((int,nonviscDampIx,-1,(Attr::hidden|Attr::noSave),"Index of the energy dissipated using the non-viscous damping (:yref:`damping<ImpulseDynamic.damping>`)."))
		((bool,kinSplit,false,,"Whether to separately track translational and rotational kinetic energy."))
		((int,kinEnergyIx,-1,(Attr::hidden|Attr::noSave),"Index for kinetic energy in scene->energies."))
		((int,kinEnergyTransIx,-1,(Attr::hidden|Attr::noSave),"Index for translational kinetic energy in scene->energies."))
		((int,kinEnergyRotIx,-1,(Attr::hidden|Attr::noSave),"Index for rotational kinetic energy in scene->energies."))
		((int,mask,-1,,"If mask defined and the bitwise AND between mask and body`s groupMask gives 0, the body will not move/rotate. Velocities and accelerations will be calculated not paying attention to this parameter."))
		,
		/*ctor*/
			densityScaling=false;
			#ifdef SUDODEM_OPENMP
				threadMaxVelocitySq.resize(omp_get_max_threads()); syncEnsured=false;
			#endif
		,/*py*/
		.add_property("densityScaling",&ImpulseDynamic::get_densityScaling,&ImpulseDynamic::set_densityScaling,"if True, then density scaling [Pfc3dManual30]_ will be applied in order to have a critical timestep equal to :yref:`GlobalStiffnessTimeStepper::targetDt` for all bodies. This option makes the simulation unrealistic from a dynamic point of view, but may speedup quasistatic simulations. In rare situations, it could be useful to not set the scalling factor automatically for each body (which the time-stepper does). In such case revert :yref:`GlobalStiffnessTimeStepper.densityScaling` to False.")
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(ImpulseDynamic);
#endif
